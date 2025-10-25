module ScriptUtils

export load_ucr_dataset, prepare_anomaly_dataset, summarize_dataset, default_time_index

using Random
using Statistics: mean
using StatsBase: countmap

struct LoadedDataset
    series::Vector{Matrix{Float64}}
    labels::Vector{String}
end

struct PreparedDataset
    series::Vector{Matrix{Float64}}
    binary_labels::Vector{Int}
    original_labels::Vector{String}
    normal_label::String
    anomaly_labels::Vector{String}
    index_map::Dict{String,Int}
end

function read_ts_file(path::AbstractString)
    lines = readlines(path)
    isempty(lines) && error("$(path) appears to be empty")
    idx = findfirst(l -> startswith(lowercase(strip(l)), "@data"), lines)
    idx === nothing && error("@data section not found in $(path)")
    data_lines = [strip(l) for l in lines[idx+1:end] if !isempty(strip(l))]
    series = Vector{Matrix{Float64}}()
    labels = String[]
    for (ln, line) in enumerate(data_lines)
        parts = split(line, ":")
        length(parts) >= 2 || error("Invalid line $(ln) in $(path): expected at least one dimension and a label")
        dims = parts[1:end-1]
        label = strip(parts[end])
        dim_vals = [parse.(Float64, filter(!isempty, split(strip(dim), ","))) for dim in dims]
        lengths = length.(dim_vals)
        all(lengths .== lengths[1]) || error("Inconsistent dimension lengths in line $(ln) of $(path)")
        n_time = lengths[1]
        n_dim = length(dim_vals)
        mat = Array{Float64}(undef, n_time, n_dim)
        for j in 1:n_dim
            mat[:, j] = dim_vals[j]
        end
        push!(series, mat)
        push!(labels, label)
    end
    LoadedDataset(series, labels)
end

function load_ucr_dataset(dataset::AbstractString, data_dir::AbstractString)
    base = joinpath(data_dir, dataset)
    train_path = base * "_TRAIN.ts"
    test_path = base * "_TEST.ts"
    isfile(train_path) || error("Training file not found: " * train_path)
    train = read_ts_file(train_path)
    test = isfile(test_path) ? read_ts_file(test_path) : LoadedDataset(Matrix{Float64}[], String[])
    LoadedDataset(vcat(train.series, test.series), vcat(train.labels, test.labels))
end

function majority_label_id(counts::Dict{Int,Int})
    max_lab = first(keys(counts))
    max_count = counts[max_lab]
    for (lab, cnt) in counts
        if cnt > max_count
            max_count = cnt
            max_lab = lab
        end
    end
    max_lab
end

function prepare_anomaly_dataset(data::LoadedDataset; anomaly_ratio::Float64 = 0.10, rng::AbstractRNG = Random.default_rng())
    length(data.series) == length(data.labels) || error("Series and labels length mismatch")
    isempty(data.series) && error("Dataset contains no series")
    label_map = Dict{String,Int}()
    mapped = Vector{Int}(undef, length(data.labels))
    next_id = 1
    for (i, lab) in enumerate(data.labels)
        lab = strip(lab)
        if !haskey(label_map, lab)
            label_map[lab] = next_id
            next_id += 1
        end
        mapped[i] = label_map[lab]
    end
    counts = countmap(mapped)
    maj_id = majority_label_id(counts)
    maj_label = ""
    for (lab, idx) in label_map
        if idx == maj_id
            maj_label = lab
            break
        end
    end
    maj_label == "" && error("Failed to identify majority class label")
    binary = [mapped[i] == maj_id ? 0 : 1 for i in eachindex(mapped)]
    normal_idx = findall(==(0), binary)
    anomaly_idx = findall(==(1), binary)
    isempty(normal_idx) && error("No majority class found in dataset")
    selected_idx = copy(normal_idx)
    if !isempty(anomaly_idx)
        n_normal = length(normal_idx)
        ratio = clamp(anomaly_ratio, 0.0, 0.99)
        n_anom_needed = ratio <= 0 ? 0 : max(1, round(Int, ratio * n_normal / max(1 - ratio, eps(Float64))))
        n_anom = min(n_anom_needed, length(anomaly_idx))
        if n_anom > 0
            shuffled_anoms = Random.shuffle(rng, anomaly_idx)
            append!(selected_idx, shuffled_anoms[1:n_anom])
        end
    else
        @warn "No anomaly samples detected; returning only normal class"
    end
    shuffled = Random.shuffle(rng, selected_idx)
    sel_series = [data.series[i] for i in shuffled]
    sel_labels = [binary[i] for i in shuffled]
    anomaly_labels = String[]
    for (lab, idx) in label_map
        idx == maj_id && continue
        push!(anomaly_labels, lab)
    end
    PreparedDataset(sel_series, sel_labels, data.labels, maj_label, anomaly_labels, label_map)
end

function summarize_dataset(prepped::PreparedDataset)
    n = length(prepped.series)
    n == 0 && return "(empty dataset)"
    n_time = size(prepped.series[1], 1)
    n_dim = size(prepped.series[1], 2)
    anomaly_pct = mean(prepped.binary_labels) * 100
    "Curves: $(n) | Length: $(n_time) | Dimensions: $(n_dim) | Anomaly %: $(round(anomaly_pct; digits=2))"
end

function default_time_index(series::Vector{Matrix{Float64}})
    isempty(series) && error("Cannot build time index for empty dataset")
    n_time = size(series[1], 1)
    collect(1:n_time)
end

end # module
