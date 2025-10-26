#!/usr/bin/env julia

# =============================================================================
# CONFIGURATION - Modify these settings as needed
# =============================================================================
# - Number of parallel workers: Set in setup_parallel_workers() function
# - MCMC iterations: Set in EXPERIMENT_CONFIG[:n_iter]
# - Parameter grid: Modify PARAMETER_GRID dictionary
# - Experiment settings: Modify EXPERIMENT_CONFIG dictionary
# =============================================================================

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))
Pkg.instantiate()

using Random
using Statistics: mean, std
using StatsBase: countmap
using TOML
using Printf
using WICMAD
using Interpolations
using Plots
using Distributed
using Base.Threads
using Wavelets
using Wavelets: WT

include(joinpath(@__DIR__, "..", "common.jl"))
using .ScriptUtils

# Make packages available on all workers (will be executed after workers are created)
# Note: These will be executed when workers are created

# Import only the specific functions we need from run_arrowhead.jl
# We'll copy the functions instead of including the whole file to avoid main() execution

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

# Function to map clusters to binary classification
function map_clusters_to_binary(clusters::Vector{Int}, truth::Vector{Int})
    normal_idx = findall(==(0), truth)
    counts = isempty(normal_idx) ? countmap(clusters) : countmap(clusters[normal_idx])
    isempty(counts) && error("Unable to determine a normal cluster from assignments")
    sorted_keys = sort(collect(keys(counts)); by = k -> -counts[k])
    normal_cluster = sorted_keys[1]
    binary = [cl == normal_cluster ? 0 : 1 for cl in clusters]
    normal_cluster, binary
end

# Function to create confusion matrix
function confusion_matrix(truth::Vector{Int}, preds::Vector{Int})
    cm = zeros(Int, 2, 2)
    for (t, p) in zip(truth, preds)
        cm[t + 1, p + 1] += 1
    end
    cm
end

# Function to select revealed indices for normal observations
# Based on the original design: reveal normal observations to help identify the normal cluster
function select_revealed_indices_normals(labels::Vector{Int}, reveal_ratio_normal::Float64, rng)
    reveal_ratio_normal <= 0 && return Int[]
    
    # Find normal observation indices (label == 0)
    normal_idx = findall(==(0), labels)
    isempty(normal_idx) && return Int[]
    
    reveal_ratio_normal = clamp(reveal_ratio_normal, 0.0, 1.0)
    n_reveal = max(1, round(Int, reveal_ratio_normal * length(normal_idx)))
    n_reveal = min(n_reveal, length(normal_idx))
    shuffle = Random.shuffle(rng, normal_idx)
    return shuffle[1:n_reveal]
end

# =============================================================================
# PARAMETER GRID CONFIGURATION
# =============================================================================
# Modify these arrays to define your factorial experiment grid

# Comprehensive wavelet families based on Wavelets.jl documentation
# Class Type: Namebase, Supertype, Numbers
const WAVELET_FAMILIES = [
    # Haar
    "haar",
    
    # Coiflet (2:2:8)
    "coif2", "coif4", "coif6", "coif8",
    
    # Daubechies (1:Inf) - using common ones
    "db1", "db2", "db3", "db4", "db5", "db6", "db7", "db8", "db9", "db10",
    
    # Symlet (4:10)
    "sym4", "sym5", "sym6", "sym7", "sym8", "sym9", "sym10",
    
    # Battle (2:2:6)
    "batt2", "batt4", "batt6",
    
    # Beylkin
    "beyl",
    
    # Vaidyanathan
    "vaid",
    
    # CDF (9,7)
    "cdf97"
]

# Model modes based on mean_intercept
# Note: WICMAD automatically sets add_bias_variants = !mean_intercept
# So we can only test 2 combinations, not 4 independently
const MODEL_MODES = [
    (mean_intercept=false, description="No mean intercept, with kernel bias variants"),
    (mean_intercept=true, description="With mean intercept, no kernel bias variants")
]

# Parameter levels for the factorial experiment
const PARAMETER_GRID = Dict(
    :alpha_prior => [(1.0, 1.0), (5.0, 1.0), (10.0, 1.0)],
    :kappa_pi => [0.1, 0.3, 0.6],
    :c2 => [0.5, 1.0, 2.0],
    :tau_pi => [10, 25, 40],
    :reveal_ratio_normal => [0.0, 0.10, 0.20],    # Fraction of normal observations revealed
    :warmup_iters => [100, 300, 500],
    :unpin => [true, false],
    :wf => WAVELET_FAMILIES,  # Now includes all wavelet families
    :model_mode => MODEL_MODES  # Different model configurations
)

# Experiment configuration
const EXPERIMENT_CONFIG = Dict(
    :n_iter => 10,            # MCMC iterations (set to 10 for quick testing)
    :burn => 5,               # Burn-in iterations
    :thin => 1,              # Thinning
    :anomaly_ratio => 0.10,  # Ratio of anomalies in dataset
    :target_length => 32,    # Target series length after interpolation
    :diagnostics => false,   # Enable/disable diagnostics for speed
    :K_init => 10           # Fixed initial number of clusters
)

# =============================================================================
# WAVELET VALIDATION
# =============================================================================

function validate_wavelets()
    """Check which wavelets are actually available in the current Wavelets.jl installation"""
    println("Validating available wavelets...")
    available_wavelets = String[]
    unavailable_wavelets = String[]
    
    for wf in WAVELET_FAMILIES
        try
            # Try to create the wavelet to see if it's available
            sym = Symbol(lowercase(wf))
            if isdefined(Wavelets.WT, sym)
                push!(available_wavelets, wf)
                println("  ✓ $wf: Available")
            else
                push!(unavailable_wavelets, wf)
                println("  ✗ $wf: Not available")
            end
        catch e
            push!(unavailable_wavelets, wf)
            println("  ✗ $wf: Error - $e")
        end
    end
    
    println("\nWavelet Validation Summary:")
    println("  Available: $(length(available_wavelets)) wavelets")
    println("  Unavailable: $(length(unavailable_wavelets)) wavelets")
    
    if !isempty(unavailable_wavelets)
        println("  Unavailable wavelets: $(join(unavailable_wavelets, ", "))")
    end
    
    return available_wavelets, unavailable_wavelets
end

# =============================================================================
# PARALLEL WORKER SETUP
# =============================================================================

function setup_parallel_workers()
    # Get number of available CPU cores
    n_cores = Sys.CPU_THREADS
    println("Available CPU cores: $n_cores")
    
    # Auto-detect number of workers (use all cores except 1 for the main process)
    n_workers = max(1, n_cores - 1)
    
    if nprocs() == 1
        addprocs(n_workers)
        println("Added $n_workers workers for parallel processing")
    else
        println("Workers already available: $(nprocs() - 1)")
    end
    
    # Print worker information
    println("Total processes: $(nprocs())")
    println("Workers: $(workers())")
    
    @everywhere begin
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

        # Function to map clusters to binary classification
        function map_clusters_to_binary(clusters::Vector{Int}, truth::Vector{Int})
            normal_idx = findall(==(0), truth)
            counts = isempty(normal_idx) ? countmap(clusters) : countmap(clusters[normal_idx])
            isempty(counts) && error("Unable to determine a normal cluster from assignments")
            sorted_keys = sort(collect(keys(counts)); by = k -> -counts[k])
            normal_cluster = sorted_keys[1]
            binary = [cl == normal_cluster ? 0 : 1 for cl in clusters]
            normal_cluster, binary
        end

        # Function to create confusion matrix
        function confusion_matrix(truth::Vector{Int}, preds::Vector{Int})
            cm = zeros(Int, 2, 2)
            for (t, p) in zip(truth, preds)
                cm[t + 1, p + 1] += 1
            end
            cm
        end

        # Function to select revealed indices for normal observations
        function select_revealed_indices_normals(labels::Vector{Int}, reveal_ratio_normal::Float64, rng)
            reveal_ratio_normal <= 0 && return Int[]
            
            # Find normal observation indices (label == 0)
            normal_idx = findall(==(0), labels)
            isempty(normal_idx) && return Int[]
            
            reveal_ratio_normal = clamp(reveal_ratio_normal, 0.0, 1.0)
            n_reveal = max(1, round(Int, reveal_ratio_normal * length(normal_idx)))
            n_reveal = min(n_reveal, length(normal_idx))
            shuffle = Random.shuffle(rng, normal_idx)
            return shuffle[1:n_reveal]
        end

        # Function to run a single experiment (will be executed on workers)
        function run_single_experiment(exp_id::Int, params::Dict, config::Dict, data_dir::String, results_file::String)
            try
                println("="^60)
                println("Worker $(myid()): Starting Experiment $exp_id")
                println("Parameters:")
                for (key, value) in params
                    println("  $key: $value")
                end
                println("="^60)
                
                # Set up random seed for reproducibility
                rng = MersenneTwister(exp_id)
                println("Worker $(myid()): Experiment $exp_id - Random seed set to $exp_id")
                
                # Load and prepare data
                println("Worker $(myid()): Experiment $exp_id - Loading ArrowHead dataset...")
                raw = ScriptUtils.load_ucr_dataset("ArrowHead", data_dir)
                prepped = ScriptUtils.prepare_anomaly_dataset(raw; 
                    anomaly_ratio = config[:anomaly_ratio], 
                    rng = rng)
                println("Worker $(myid()): Experiment $exp_id - Dataset prepared: $(length(prepped.series)) series")
                
                # Interpolate to target length
                t = ScriptUtils.default_time_index(prepped.series)
                original_length = length(t)
                target_length = config[:target_length]
                
                if original_length != target_length
                    println("Worker $(myid()): Experiment $exp_id - Interpolating from $original_length to $target_length points")
                    t_old = collect(1.0:original_length)
                    t_new = collect(range(1.0, original_length, length=target_length))
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
                else
                    println("Worker $(myid()): Experiment $exp_id - No interpolation needed (length = $original_length)")
                end
                
                # Select revealed indices for normal observations only
                reveal_idx = select_revealed_indices_normals(prepped.binary_labels, params[:reveal_ratio_normal], rng)
                println("Worker $(myid()): Experiment $exp_id - Revealed $(length(reveal_idx)) normal observations")
                
                # Extract model mode parameters
                model_mode = params[:model_mode]
                mean_intercept = model_mode.mean_intercept
                
                # Note: WICMAD automatically sets add_bias_variants = !mean_intercept
                # So we can only test 2 combinations, not 4 independently
                
                println("Worker $(myid()): Experiment $exp_id - Starting WICMAD with $(config[:n_iter]) iterations...")
                println("Worker $(myid()): Experiment $exp_id - Full config dictionary:")
                for (key, value) in config
                    println("  $key: $value")
                end
                println("Worker $(myid()): Experiment $exp_id - WICMAD parameters:")
                println("  n_iter: $(config[:n_iter])")
                println("  burn: $(config[:burn])")
                println("  thin: $(config[:thin])")
                println("  wf: $(params[:wf])")
                println("  mean_intercept: $mean_intercept")
                start_time = time()
                
                # Run WICMAD with custom parameters
                result = wicmad(prepped.series, t;
                    n_iter = config[:n_iter],
                    burn = config[:burn],
                    thin = config[:thin],
                    alpha_prior = params[:alpha_prior],
                    wf = params[:wf],
                    kappa_pi = params[:kappa_pi],
                    c2 = params[:c2],
                    tau_pi = Int(params[:tau_pi]),
                    warmup_iters = params[:warmup_iters],
                    unpin = params[:unpin],
                    K_init = config[:K_init],
                    revealed_idx = reveal_idx,
                    diagnostics = config[:diagnostics],
                    mean_intercept = mean_intercept
                )
                
                wicmad_time = time() - start_time
                println("Worker $(myid()): Experiment $exp_id - WICMAD completed in $(round(wicmad_time, digits=2)) seconds")
                
                # Process results
                println("Worker $(myid()): Experiment $exp_id - Processing results...")
                Z = result.Z
                samples = size(Z, 1)
                clusters = samples > 0 ? vec(Z[end, :]) : ones(Int, length(prepped.series))
                truth = prepped.binary_labels
                normal_cluster, binary_preds = map_clusters_to_binary(clusters, truth)
                
                # Calculate metrics
                ari = WICMAD.adj_rand_index(binary_preds, truth)
                cm = confusion_matrix(truth, binary_preds)
                
                # Calculate F1 score for anomaly class (class 1)
                tp = sum((truth .== 1) .& (binary_preds .== 1))
                fp = sum((truth .== 0) .& (binary_preds .== 1))
                fn = sum((truth .== 1) .& (binary_preds .== 0))
                precision = (tp + fp) > 0 ? tp / (tp + fp) : 0.0
                recall = (tp + fn) > 0 ? tp / (tp + fn) : 0.0
                f1_score = (precision + recall) > 0 ? 2 * precision * recall / (precision + recall) : 0.0
                
                # Calculate overall accuracy
                accuracy = sum(truth .== binary_preds) / length(truth)
                
                println("Worker $(myid()): Experiment $exp_id - Results calculated:")
                println("  F1 Score: $(round(f1_score, digits=4))")
                println("  Accuracy: $(round(accuracy, digits=4))")
                println("  ARI: $(round(ari, digits=4))")
                println("  Clusters found: $(length(unique(clusters)))")
                println("  Normal cluster: $normal_cluster")
                
                # Record results
                result_data = Dict(
                    "experiment_id" => exp_id,
                    "worker_id" => myid(),
                    "alpha_prior" => params[:alpha_prior],
                    "kappa_pi" => params[:kappa_pi],
                    "c2" => params[:c2],
                    "tau_pi" => params[:tau_pi],
                    "reveal_ratio_normal" => params[:reveal_ratio_normal],
                    "warmup_iters" => params[:warmup_iters],
                    "unpin" => params[:unpin],
                    "K_init" => config[:K_init],
                    "wf" => params[:wf],
                    "mean_intercept" => mean_intercept,
                    "model_description" => model_mode.description,
                    "f1_score" => f1_score,
                    "precision" => precision,
                    "recall" => recall,
                    "accuracy" => accuracy,
                    "ari" => ari,
                    "n_clusters" => length(unique(clusters)),
                    "normal_cluster" => normal_cluster,
                    "tp" => tp,
                    "fp" => fp,
                    "fn" => fn,
                    "tn" => sum((truth .== 0) .& (binary_preds .== 0))
                )
                
                # Write result to file (thread-safe)
                println("Worker $(myid()): Experiment $exp_id - Writing results to file...")
                open(results_file, "a") do io
                    TOML.print(io, Dict("experiment_$exp_id" => result_data))
                    println(io)  # Add newline
                end
                
                println("="^60)
                println("Worker $(myid()): Experiment $exp_id COMPLETED SUCCESSFULLY!")
                println("Final Results: F1 = $(round(f1_score, digits=4)), Accuracy = $(round(accuracy, digits=4))")
                println("="^60)
                return result_data
                
            catch e
                println("Worker $(myid()): Experiment $exp_id failed: $e")
                error_data = Dict(
                    "experiment_id" => exp_id,
                    "worker_id" => myid(),
                    "error" => string(e),
                    "f1_score" => 0.0,
                    "accuracy" => 0.0
                )
                
                # Write error to file
                open(results_file, "a") do io
                    TOML.print(io, Dict("experiment_$exp_id" => error_data))
                    println(io)
                end
                
                return error_data
            end
        end
    end
    
    return n_workers
end

# =============================================================================
# EXPERIMENT FUNCTIONS
# =============================================================================

# Function to create parameter combinations from the grid
function create_parameter_combinations(grid::Dict)
    parameter_combinations = []
    exp_id = 1
    
    # Get all parameter names and their values
    param_names = collect(keys(grid))
    param_values = collect(values(grid))
    
    # Create all combinations using recursive function
    function generate_combinations(current_params::Dict, remaining_params::Vector)
        if isempty(remaining_params)
            push!(parameter_combinations, (exp_id, copy(current_params)))
            return exp_id + 1
        else
            param_name, param_values = popfirst!(remaining_params)
            for value in param_values
                current_params[param_name] = value
                exp_id = generate_combinations(current_params, copy(remaining_params))
            end
        end
        return exp_id
    end
    
    # Convert grid to vector of pairs for easier manipulation
    param_pairs = [name => values for (name, values) in grid]
    generate_combinations(Dict(), param_pairs)
    
    return parameter_combinations
end

# Function to analyze results
function analyze_results(results_file::String)
    println("\n" * "="^80)
    println("PARALLEL FACTORIAL EXPERIMENT RESULTS")
    println("="^80)
    
    # Read all results
    if !isfile(results_file)
        println("Results file not found: $results_file")
        return nothing
    end
    
    results = TOML.parsefile(results_file)
    
    # Extract F1 scores and parameters
    f1_scores = Float64[]
    experiment_data = []
    worker_stats = Dict{Int, Int}()
    
    for (key, data) in results
        if startswith(key, "experiment_")
            f1 = get(data, "f1_score", 0.0)
            push!(f1_scores, f1)
            push!(experiment_data, data)
            
            # Track worker usage
            worker_id = get(data, "worker_id", 1)
            worker_stats[worker_id] = get(worker_stats, worker_id, 0) + 1
        end
    end
    
    if isempty(f1_scores)
        println("No valid results found!")
        return nothing
    end
    
    # Find best experiment
    best_idx = argmax(f1_scores)
    best_experiment = experiment_data[best_idx]
    best_f1 = f1_scores[best_idx]
    
    println("Total experiments: $(length(f1_scores))")
    println("Best F1 Score: $(round(best_f1, digits=4))")
    println("Mean F1 Score: $(round(mean(f1_scores), digits=4))")
    println("Std F1 Score: $(round(std(f1_scores), digits=4))")
    
    # Worker statistics
    println("\nWorker Usage Statistics:")
    println("-"^40)
    for (worker_id, count) in sort(collect(worker_stats))
        println("  Worker $worker_id: $count experiments")
    end
    
    println("\nBest Parameterization:")
    println("-"^40)
    for (param, value) in best_experiment
        if param != "experiment_id" && param != "worker_id" && param != "f1_score" && 
           param != "precision" && param != "recall" && param != "accuracy" && 
           param != "ari" && param != "n_clusters" && param != "normal_cluster" && 
           param != "tp" && param != "fp" && param != "fn" && param != "tn" && 
           param != "error" && param != "model_description"
            println("  $param: $value")
        end
    end
    
    # Show model description separately
    if haskey(best_experiment, "model_description")
        println("  Model Mode: $(best_experiment["model_description"])")
    end
    
    println("\nBest Experiment Performance:")
    println("-"^40)
    println("  F1 Score: $(round(best_experiment["f1_score"], digits=4))")
    println("  Precision: $(round(best_experiment["precision"], digits=4))")
    println("  Recall: $(round(best_experiment["recall"], digits=4))")
    println("  Accuracy: $(round(best_experiment["accuracy"], digits=4))")
    println("  ARI: $(round(best_experiment["ari"], digits=4))")
    
    # Top 5 experiments
    sorted_indices = sortperm(f1_scores, rev=true)
    println("\nTop 5 Experiments:")
    println("-"^40)
    for i in 1:min(5, length(sorted_indices))
        idx = sorted_indices[i]
        exp_data = experiment_data[idx]
        println("  #$i: F1=$(round(exp_data["f1_score"], digits=4)), " *
                "ExpID=$(exp_data["experiment_id"]), " *
                "Worker=$(exp_data["worker_id"]), " *
                "wf=$(exp_data["wf"]), " *
                "K_init=$(exp_data["K_init"])")
    end
    
    # Model mode analysis
    println("\nModel Mode Performance Analysis:")
    println("-"^40)
    model_stats = Dict{String, Vector{Float64}}()
    for exp_data in experiment_data
        model_desc = exp_data["model_description"]
        f1 = exp_data["f1_score"]
        if !haskey(model_stats, model_desc)
            model_stats[model_desc] = Float64[]
        end
        push!(model_stats[model_desc], f1)
    end
    
    println("Model Mode Performance:")
    for (model_desc, scores) in model_stats
        mean_f1 = mean(scores)
        std_f1 = std(scores)
        count = length(scores)
        println("  $model_desc:")
        println("    Mean F1: $(round(mean_f1, digits=4)), Std: $(round(std_f1, digits=4)), Count: $count")
    end
    
    # Wavelet-specific analysis
    println("\nWavelet Performance Analysis:")
    println("-"^40)
    wavelet_stats = Dict{String, Vector{Float64}}()
    for exp_data in experiment_data
        wf = exp_data["wf"]
        f1 = exp_data["f1_score"]
        if !haskey(wavelet_stats, wf)
            wavelet_stats[wf] = Float64[]
        end
        push!(wavelet_stats[wf], f1)
    end
    
    # Sort wavelets by mean F1 score
    wavelet_means = [(wf, mean(scores)) for (wf, scores) in wavelet_stats]
    sort!(wavelet_means, by=x->x[2], rev=true)
    
    println("Top 10 Wavelets by Mean F1 Score:")
    for i in 1:min(10, length(wavelet_means))
        wf, mean_f1 = wavelet_means[i]
        scores = wavelet_stats[wf]
        println("  #$i: $wf - Mean F1: $(round(mean_f1, digits=4)), " *
                "Std: $(round(std(scores), digits=4)), " *
                "Count: $(length(scores))")
    end
    
    return best_experiment
end

# =============================================================================
# MAIN EXECUTION
# =============================================================================

function setup_workers_and_packages()
    # Setup parallel workers
    n_workers = setup_parallel_workers()
    return n_workers
end

function main()
    # Open output file for all print statements
    output_file = open("parallel_factorial_output.txt", "w")
    
    # Redirect stdout to file
    redirect_stdout(output_file) do
        println("PARALLEL FACTORIAL EXPERIMENT")
        println("="^50)
    
    # Validate wavelets first
    available_wavelets, unavailable_wavelets = validate_wavelets()
    
    # Update parameter grid to only use available wavelets
    if !isempty(unavailable_wavelets)
        println("\nUpdating parameter grid to use only available wavelets...")
        PARAMETER_GRID[:wf] = available_wavelets
    end
    
    # Setup parallel workers
    n_workers = setup_workers_and_packages()
    
    # Load packages on all workers using eval
    for w in workers()
        remotecall_wait(w) do
            eval(:(using Random))
            eval(:(using Statistics: mean, std))
            eval(:(using StatsBase: countmap))
            eval(:(using TOML))
            eval(:(using Printf))
            eval(:(using WICMAD))
            eval(:(using Interpolations))
            eval(:(using Plots))
            eval(:(using Wavelets))
            eval(:(using Wavelets: WT))
            
            eval(:(include(joinpath(@__DIR__, "..", "common.jl"))))
            eval(:(using .ScriptUtils))
        end
    end
    
    # Setup paths
    data_dir = joinpath(@__DIR__, "..", "..", "..", "data", "ArrowHead")
    results_file = joinpath(@__DIR__, "parallel_factorial_results.toml")
    
    # Clear results file
    if isfile(results_file)
        rm(results_file)
        println("Cleared previous results file")
    end
    
    # Create parameter combinations
    parameter_combinations = create_parameter_combinations(PARAMETER_GRID)
    total_experiments = length(parameter_combinations)
    
    println("\nExperiment Configuration:")
    println("-"^40)
    println("Total experiments: $total_experiments")
    println("Results file: $results_file")
    println("Data directory: $data_dir")
    println("Workers: $n_workers")
    
    # Print parameter grid summary
    println("\nParameter Grid:")
    println("-"^40)
    for (param, values) in PARAMETER_GRID
        println("  $param: $(length(values)) levels $(values)")
    end
    
    # Check if data directory exists
    if !isdir(data_dir)
        println("ERROR: Data directory not found: $data_dir")
        return
    end
    
    # Run experiments in parallel
    println("\nStarting parallel experiments...")
    println("Progress will be shown by each worker as experiments complete.")
    start_time = time()
    
    # Use pmap for parallel execution across workers
    println("Launching $(total_experiments) experiments across $(n_workers) workers...")
    results = pmap(parameter_combinations) do (exp_id, params)
        run_single_experiment(exp_id, params, EXPERIMENT_CONFIG, data_dir, results_file)
    end
    
    end_time = time()
    elapsed_time = end_time - start_time
    
    println("\nAll experiments completed!")
    println("Total time: $(round(elapsed_time, digits=2)) seconds")
    println("Average time per experiment: $(round(elapsed_time / total_experiments, digits=2)) seconds")
    println("Speedup with $n_workers workers: $(round(total_experiments * elapsed_time / total_experiments / n_workers, digits=2))x theoretical")
    
    # Analyze results
    analyze_results(results_file)
    
    # Clean up workers
    rmprocs(workers())
    println("\nWorkers removed. Experiment complete.")
    end  # End redirect_stdout block
    
    # Close output file
    close(output_file)
    println("All output saved to parallel_factorial_output.txt")
end

# Run the experiment immediately
main()
