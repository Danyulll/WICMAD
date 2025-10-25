#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using Random
using Statistics: mean
using TOML
using WICMAD

include(joinpath(@__DIR__, "common.jl"))
using .ScriptUtils

struct Config
    dataset::String
    data_dir::String
    anomaly_ratio::Float64
    reveal_ratio::Float64
    seed::Int
    n_iter::Int
    burn::Int
    thin::Int
    diagnostics::Bool
    mean_intercept::Bool
    output_path::Union{Nothing,String}
end

function usage()
    println("Usage: julia --project=julia scripts/run_wicmad.jl --dataset DATASET --data-dir PATH [options]\n")
    println("Options:")
    println("  --dataset NAME           Dataset prefix (e.g. ArrowHead)")
    println("  --data-dir PATH          Directory containing NAME_TRAIN.ts and NAME_TEST.ts")
    println("  --anomaly-ratio R        Target anomaly ratio in (0,1); default 0.10")
    println("  --reveal-ratio R         Fraction of anomalies to reveal to WICMAD; default 0.15")
    println("  --n-iter N               Total MCMC iterations (default 6000)")
    println("  --burn N                 Burn-in iterations (default 3000)")
    println("  --thin N                 Thinning interval (default 5)")
    println("  --seed N                 RNG seed (default 1)")
    println("  --mean-intercept         Enable mean-intercept sampling")
    println("  --no-diagnostics         Disable diagnostic collection")
    println("  --output PATH            Optional TOML file to record summary metrics")
    println("  --help                   Show this message")
end

function parse_args(args::Vector{String})
    dataset = nothing
    data_dir = pwd()
    anomaly_ratio = 0.10
    reveal_ratio = 0.15
    seed = 1
    n_iter = 6000
    burn = 3000
    thin = 5
    diagnostics = true
    mean_intercept = false
    output_path = nothing
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg in ("--dataset", "-d")
            i += 1
            i > length(args) && error("--dataset requires a value")
            dataset = args[i]
        elseif arg in ("--data-dir", "-D")
            i += 1
            i > length(args) && error("--data-dir requires a value")
            data_dir = args[i]
        elseif arg == "--anomaly-ratio"
            i += 1
            i > length(args) && error("--anomaly-ratio requires a value")
            anomaly_ratio = parse(Float64, args[i])
        elseif arg == "--reveal-ratio"
            i += 1
            i > length(args) && error("--reveal-ratio requires a value")
            reveal_ratio = parse(Float64, args[i])
        elseif arg == "--seed"
            i += 1
            i > length(args) && error("--seed requires a value")
            seed = parse(Int, args[i])
        elseif arg == "--n-iter"
            i += 1
            i > length(args) && error("--n-iter requires a value")
            n_iter = parse(Int, args[i])
        elseif arg == "--burn"
            i += 1
            i > length(args) && error("--burn requires a value")
            burn = parse(Int, args[i])
        elseif arg == "--thin"
            i += 1
            i > length(args) && error("--thin requires a value")
            thin = parse(Int, args[i])
        elseif arg == "--mean-intercept"
            mean_intercept = true
        elseif arg == "--no-diagnostics"
            diagnostics = false
        elseif arg == "--output"
            i += 1
            i > length(args) && error("--output requires a value")
            output_path = args[i]
        elseif arg in ("--help", "-h")
            usage()
            exit(0)
        else
            error("Unknown argument: $(arg)")
        end
        i += 1
    end
    dataset === nothing && error("--dataset must be provided")
    Config(dataset, data_dir, anomaly_ratio, reveal_ratio, seed, n_iter, burn, thin, diagnostics, mean_intercept, output_path)
end

function select_revealed_indices(labels::Vector{Int}, reveal_ratio::Float64, rng)
    reveal_ratio <= 0 && return Int[]
    idx = findall(==(1), labels)
    isempty(idx) && return Int[]
    reveal_ratio = clamp(reveal_ratio, 0.0, 1.0)
    n_reveal = max(1, round(Int, reveal_ratio * length(idx)))
    n_reveal = min(n_reveal, length(idx))
    shuffle = Random.shuffle(rng, idx)
    shuffle[1:n_reveal]
end

function main(args)
    cfg = parse_args(args)
    rng = MersenneTwister(cfg.seed)
    println("Loading dataset $(cfg.dataset) from $(cfg.data_dir)...")
    raw = load_ucr_dataset(cfg.dataset, cfg.data_dir)
    prepped = prepare_anomaly_dataset(raw; anomaly_ratio = cfg.anomaly_ratio, rng = rng)
    println("Prepared dataset summary: " * summarize_dataset(prepped))
    println("Normal class: " * prepped.normal_label)
    if !isempty(prepped.anomaly_labels)
        println("Anomaly classes: " * join(prepped.anomaly_labels, ", "))
    else
        println("Anomaly classes: (none)")
    end
    t = default_time_index(prepped.series)
    reveal_idx = select_revealed_indices(prepped.binary_labels, cfg.reveal_ratio, rng)
    if !isempty(reveal_idx)
        println("Revealing $(length(reveal_idx)) anomalies to the sampler")
    end
    result = wicmad(prepped.series, t; n_iter = cfg.n_iter, burn = cfg.burn, thin = cfg.thin,
        diagnostics = cfg.diagnostics, revealed_idx = reveal_idx, mean_intercept = cfg.mean_intercept)
    Z = result.Z
    est_labels = size(Z, 1) > 0 ? vec(Z[end, :]) : ones(Int, length(prepped.series))
    truth = [prepped.binary_labels[i] == 0 ? 1 : 2 for i in eachindex(prepped.binary_labels)]
    ari = WICMAD.adj_rand_index(est_labels, truth)
    uniq_clusters = length(unique(est_labels))
    println("\nWICMAD run complete:")
    println("  Samples saved: " * string(size(Z, 1)))
    println("  Unique clusters (last sample): " * string(uniq_clusters))
    println("  Adjusted Rand Index vs binary truth: " * string(round(ari; digits = 4)))
    if cfg.output_path !== nothing
        summary = Dict(
            "dataset" => cfg.dataset,
            "data_dir" => cfg.data_dir,
            "curves" => length(prepped.series),
            "length" => size(prepped.series[1], 1),
            "dimensions" => size(prepped.series[1], 2),
            "anomaly_ratio" => mean(prepped.binary_labels),
            "revealed" => length(reveal_idx),
            "unique_clusters" => uniq_clusters,
            "ari" => ari,
            "n_iter" => cfg.n_iter,
            "burn" => cfg.burn,
            "thin" => cfg.thin,
            "diagnostics" => cfg.diagnostics,
            "mean_intercept" => cfg.mean_intercept,
        )
        open(cfg.output_path, "w") do io
            TOML.print(io, summary)
        end
        println("Summary written to " * cfg.output_path)
    end
end

isinteractive() || main(ARGS)
