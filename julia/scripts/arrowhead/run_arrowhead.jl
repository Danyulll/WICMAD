#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))
Pkg.instantiate()

using Random
using Statistics: mean
using StatsBase: countmap
using TOML
using Printf
using WICMAD
using Interpolations

include(joinpath(@__DIR__, "..", "common.jl"))
using .ScriptUtils

# Linear interpolation function for time series using Interpolations.jl
function interpolate_series(series::Matrix{Float64}, t_old::Vector{Float64}, t_new::Vector{Float64})
    n_dims = size(series, 2)
    n_new = length(t_new)
    new_series = zeros(Float64, n_new, n_dims)
    
    for dim in 1:n_dims
        # Create linear interpolation object
        itp = LinearInterpolation(t_old, series[:, dim], extrapolation_bc=Line())
        # Evaluate at new time points
        new_series[:, dim] = itp.(t_new)
    end
    
    return new_series
end

struct ArrowheadConfig
    data_dir::String
    anomaly_ratio::Float64
    reveal_ratio::Float64
    seed::Int
    n_iter::Int
    burn::Int
    thin::Int
    warmup::Int
    diagnostics::Bool
    mean_intercept::Bool
    metrics_path::Union{Nothing,String}
end

function usage()
    println("Usage: julia --project=julia scripts/arrowhead/run_arrowhead.jl [options]\n")
    println("Options:")
    println("  --data-dir PATH          Directory containing ArrowHead_TRAIN.ts and ArrowHead_TEST.ts")
    println("  --anomaly-ratio R        Target anomaly ratio in (0,1); default 0.10")
    println("  --reveal-ratio R         Fraction of anomalies to reveal to the sampler; default 0.15")
    println("  --seed N                 RNG seed (default 1)")
    println("  --n-iter N               Total MCMC iterations (default 100)")
    println("  --burn N                 Burn-in iterations (default 50)")
    println("  --thin N                 Thinning interval (default 1)")
    println("  --warmup N               Warmup iterations passed to WICMAD (default 500)")
    println("  --mean-intercept         Enable mean-intercept sampling")
    println("  --no-diagnostics         Disable diagnostic collection")
    println("  --metrics PATH           Optional TOML file to record dataset metrics")
    println("  --help                   Show this message")
end

function parse_args(args::Vector{String})
    # Try to find the data directory relative to the project root
    project_root = joinpath(@__DIR__, "..", "..", "..")
    default_dir = joinpath(project_root, "data", "ArrowHead")
    data_dir = default_dir
    anomaly_ratio = 0.10
    reveal_ratio = 0.15
    seed = 1
    n_iter = 5000
    burn = 2000
    thin = 1
    warmup = 500
    diagnostics = true
    mean_intercept = false
    metrics_path = nothing
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg in ("--data-dir", "-D")
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
        elseif arg == "--warmup"
            i += 1
            i > length(args) && error("--warmup requires a value")
            warmup = parse(Int, args[i])
        elseif arg == "--mean-intercept"
            mean_intercept = true
        elseif arg == "--no-diagnostics"
            diagnostics = false
        elseif arg == "--metrics"
            i += 1
            i > length(args) && error("--metrics requires a value")
            metrics_path = args[i]
        elseif arg in ("--help", "-h")
            usage()
            exit(0)
        else
            error("Unknown argument: $(arg)")
        end
        i += 1
    end
    if data_dir == default_dir && !isdir(data_dir)
        @warn "Default ArrowHead data directory not found: $data_dir"
        @warn "Please specify --data-dir PATH to point to the directory containing ArrowHead_TRAIN.ts and ArrowHead_TEST.ts"
        @warn "Example: --data-dir \"C:\\Users\\danie\\Desktop\\WICMAD\\data\\ArrowHead\""
    end
    ArrowheadConfig(data_dir, anomaly_ratio, reveal_ratio, seed, n_iter, burn, thin, warmup, diagnostics, mean_intercept, metrics_path)
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

function map_clusters_to_binary(clusters::Vector{Int}, truth::Vector{Int})
    normal_idx = findall(==(0), truth)
    counts = isempty(normal_idx) ? countmap(clusters) : countmap(clusters[normal_idx])
    isempty(counts) && error("Unable to determine a normal cluster from assignments")
    sorted_keys = sort(collect(keys(counts)); by = k -> -counts[k])
    normal_cluster = sorted_keys[1]
    binary = [cl == normal_cluster ? 0 : 1 for cl in clusters]
    normal_cluster, binary
end

safe_div(num, denom) = denom == 0 ? 0.0 : num / denom

function class_metric(truth::Vector{Int}, preds::Vector{Int}, cls::Int)
    tp = sum((truth .== cls) .& (preds .== cls))
    fp = sum((truth .!= cls) .& (preds .== cls))
    fn = sum((truth .== cls) .& (preds .!= cls))
    tn = length(truth) - tp - fp - fn
    precision = safe_div(tp, tp + fp)
    recall = safe_div(tp, tp + fn)
    specificity = safe_div(tn, tn + fp)
    f1 = safe_div(2 * precision * recall, precision + recall)
    accuracy = safe_div(tp + tn, tp + tn + fp + fn)
    Dict(
        :precision => precision,
        :recall => recall,
        :specificity => specificity,
        :f1 => f1,
        :accuracy => accuracy,
        :tp => tp,
        :tn => tn,
        :fp => fp,
        :fn => fn,
    )
end

function confusion_matrix(truth::Vector{Int}, preds::Vector{Int})
    cm = zeros(Int, 2, 2)
    for (t, p) in zip(truth, preds)
        cm[t + 1, p + 1] += 1
    end
    cm
end

function macro_metrics(metrics::Dict{Int,Dict{Symbol,Float64}})
    isempty(metrics) && return Dict{Symbol,Float64}()
    classes = Set(keys(metrics))
    classes == Set([0, 1]) || return Dict{Symbol,Float64}()
    precision = mean([metrics[c][:precision] for c in classes])
    recall = mean([metrics[c][:recall] for c in classes])
    specificity = mean([metrics[c][:specificity] for c in classes])
    f1 = mean([metrics[c][:f1] for c in classes])
    accuracy = mean([metrics[c][:accuracy] for c in classes])
    Dict(
        :macro_precision => precision,
        :macro_recall => recall,
        :macro_specificity => specificity,
        :macro_f1 => f1,
        :macro_accuracy => accuracy,
    )
end

function print_metrics(metrics_per_class::Dict{Int,Dict{Symbol,Float64}}, macro_metrics::Dict{Symbol,Float64}, cm)
    println("\nBinary classification metrics (truth rows x predicted columns):")
    println("  Confusion matrix:")
    println(@sprintf("      Truth=0 -> [%4d %4d]", cm[1, 1], cm[1, 2]))
    println(@sprintf("      Truth=1 -> [%4d %4d]", cm[2, 1], cm[2, 2]))
    for cls in sort(collect(keys(metrics_per_class)))
        m = metrics_per_class[cls]
        println("\n  Class $(cls) metrics:")
        println(@sprintf("    Precision: %.4f", m[:precision]))
        println(@sprintf("    Recall: %.4f", m[:recall]))
        println(@sprintf("    Specificity: %.4f", m[:specificity]))
        println(@sprintf("    F1 score: %.4f", m[:f1]))
        println(@sprintf("    Accuracy: %.4f", m[:accuracy]))
    end
    if !isempty(macro_metrics)
        println("\n  Macro averages:")
        println(@sprintf("    Precision: %.4f", macro_metrics[:macro_precision]))
        println(@sprintf("    Recall: %.4f", macro_metrics[:macro_recall]))
        println(@sprintf("    Specificity: %.4f", macro_metrics[:macro_specificity]))
        println(@sprintf("    F1 score: %.4f", macro_metrics[:macro_f1]))
        println(@sprintf("    Accuracy: %.4f", macro_metrics[:macro_accuracy]))
    end
end

function record_metrics(path::String, cfg::ArrowheadConfig, summary)
    dir = dirname(path)
    isdir(dir) || mkpath(dir)
    open(path, "w") do io
        TOML.print(io, summary)
    end
    println("Metrics written to " * path)
end

function main(args)
    cfg = parse_args(args)
    rng = MersenneTwister(cfg.seed)
    
    # Check if data directory exists before trying to load
    if !isdir(cfg.data_dir)
        println("ERROR: Data directory not found: $(cfg.data_dir)")
        println("Please run the script with --data-dir argument pointing to the ArrowHead data directory.")
        println("Example: julia run_arrowhead.jl --data-dir \"C:\\Users\\danie\\Desktop\\WICMAD\\data\\ArrowHead\"")
        return
    end
    
    println("Loading ArrowHead dataset from $(cfg.data_dir)...")
    raw = ScriptUtils.load_ucr_dataset("ArrowHead", cfg.data_dir)
    prepped = ScriptUtils.prepare_anomaly_dataset(raw; anomaly_ratio = cfg.anomaly_ratio, rng = rng)
    println("Prepared dataset summary: " * ScriptUtils.summarize_dataset(prepped))
    println("Normal class label: " * prepped.normal_label)
    if isempty(prepped.anomaly_labels)
        println("Anomaly class labels: (none)")
    else
        println("Anomaly class labels: " * join(prepped.anomaly_labels, ", "))
    end
    t = ScriptUtils.default_time_index(prepped.series)
    
    # Interpolate data to 32 points for WICMAD compatibility
    original_length = length(t)
    target_length = 32
    if original_length != target_length
        println("Interpolating time series from $original_length to $target_length points for WICMAD compatibility")
        
        # Create new time indices for interpolation
        t_old = collect(1.0:original_length)
        t_new = collect(range(1.0, original_length, length=target_length))
        
        # Interpolate each series
        interpolated_series = [interpolate_series(series, t_old, t_new) for series in prepped.series]
        prepped = ScriptUtils.PreparedDataset(
            interpolated_series,
            prepped.binary_labels,
            prepped.original_labels,
            prepped.normal_label,
            prepped.anomaly_labels,
            prepped.index_map
        )
        t = t_new
    end
    
    reveal_idx = select_revealed_indices(prepped.binary_labels, cfg.reveal_ratio, rng)
    if !isempty(reveal_idx)
        println("Revealing $(length(reveal_idx)) anomalies to the sampler")
    end
    println("\nRunning WICMAD with n_iter=$(cfg.n_iter), burn=$(cfg.burn), thin=$(cfg.thin)...")
    println("="^60)
    result = wicmad(prepped.series, t;
        n_iter = cfg.n_iter,
        burn = cfg.burn,
        thin = cfg.thin,
        warmup_iters = cfg.warmup,
        diagnostics = cfg.diagnostics,
        revealed_idx = reveal_idx,
        mean_intercept = cfg.mean_intercept,
        wf = "haar",  # Use Haar wavelet
    )
    println("="^60)
    println("WICMAD algorithm completed!")
    Z = result.Z
    samples = size(Z, 1)
    println("Samples retained: $(samples)")
    clusters = samples > 0 ? vec(Z[end, :]) : ones(Int, length(prepped.series))
    truth = prepped.binary_labels
    normal_cluster, binary_preds = map_clusters_to_binary(clusters, truth)
    println("Normal cluster identified as label $(normal_cluster)")
    ari = WICMAD.adj_rand_index(binary_preds, truth)
    println(@sprintf("Adjusted Rand Index (binary mapping vs truth): %.4f", ari))
    cm = confusion_matrix(truth, binary_preds)
    metrics_per_class = Dict{Int,Dict{Symbol,Float64}}()
    for cls in sort(unique(truth))
        metrics_per_class[cls] = class_metric(truth, binary_preds, cls)
    end
    macro_metrics_result = macro_metrics(metrics_per_class)
    print_metrics(metrics_per_class, macro_metrics_result, cm)
    if cfg.metrics_path !== nothing
        summary = Dict(
            "dataset" => "ArrowHead",
            "data_dir" => cfg.data_dir,
            "curves" => length(prepped.series),
            "length" => size(prepped.series[1], 1),
            "dimensions" => size(prepped.series[1], 2),
            "anomaly_ratio" => mean(prepped.binary_labels),
            "revealed" => length(reveal_idx),
            "n_iter" => cfg.n_iter,
            "burn" => cfg.burn,
            "thin" => cfg.thin,
            "warmup" => cfg.warmup,
            "mean_intercept" => cfg.mean_intercept,
            "diagnostics" => cfg.diagnostics,
            "samples_kept" => samples,
            "normal_cluster" => normal_cluster,
            "ari" => ari,
            "confusion_matrix" => [[cm[1,1], cm[1,2]], [cm[2,1], cm[2,2]]],
            "class_metrics" => Dict(string(k) => Dict(string(sym) => metrics_per_class[k][sym] for sym in keys(metrics_per_class[k])) for k in keys(metrics_per_class)),
            "macro_metrics" => Dict(string(sym) => macro_metrics_result[sym] for sym in keys(macro_metrics_result)),
        )
        record_metrics(cfg.metrics_path, cfg, summary)
    end
end

main(ARGS)
