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
using Plots
using Distributed
using SharedArrays

include(joinpath(@__DIR__, "..", "common.jl"))
using .ScriptUtils

# Import the run_arrowhead functions
include("run_arrowhead.jl")

# Function to run a single experiment
function run_single_experiment(exp_id::Int, params::Dict, data_dir::String, results_file::String)
    try
        println("Experiment $exp_id: Starting with parameters: $params")
        
        # Set up random seed for reproducibility
        rng = MersenneTwister(exp_id)
        
        # Load and prepare data
        raw = ScriptUtils.load_ucr_dataset("ArrowHead", data_dir)
        prepped = ScriptUtils.prepare_anomaly_dataset(raw; 
            anomaly_ratio = 0.10, 
            rng = rng)
        
        # Interpolate to 32 points
        t = ScriptUtils.default_time_index(prepped.series)
        original_length = length(t)
        target_length = 32
        if original_length != target_length
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
        end
        
        # Select revealed indices
        reveal_idx = select_revealed_indices(prepped.binary_labels, params[:reveal_ratio], rng)
        
        # Run WICMAD with custom parameters
        result = wicmad(prepped.series, t;
            n_iter = 100,  # Reduced for faster experiments
            burn = 50,
            thin = 1,
            alpha_prior = params[:alpha_prior],
            wf = params[:wf],
            kappa_pi = params[:kappa_pi],
            c2 = params[:c2],
            tau_pi = params[:tau_pi],
            warmup_iters = params[:warmup_iters],
            unpin = params[:unpin],
            K_init = params[:K_init],
            revealed_idx = reveal_idx,
            diagnostics = false,  # Disable for speed
            mean_intercept = false
        )
        
        # Process results
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
        
        # Record results
        result_data = Dict(
            "experiment_id" => exp_id,
            "alpha_prior" => params[:alpha_prior],
            "kappa_pi" => params[:kappa_pi],
            "c2" => params[:c2],
            "tau_pi" => params[:tau_pi],
            "reveal_ratio" => params[:reveal_ratio],
            "warmup_iters" => params[:warmup_iters],
            "unpin" => params[:unpin],
            "K_init" => params[:K_init],
            "wf" => params[:wf],
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
        open(results_file, "a") do io
            TOML.print(io, Dict("experiment_$exp_id" => result_data))
            println(io)  # Add newline
        end
        
        println("Experiment $exp_id: F1 = $(round(f1_score, digits=4)), Accuracy = $(round(accuracy, digits=4))")
        return result_data
        
    catch e
        println("Experiment $exp_id failed: $e")
        error_data = Dict(
            "experiment_id" => exp_id,
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

# Function to create parameter grid
function create_parameter_grid()
    # Define parameter levels
    alpha_prior_levels = [(1.0, 1.0), (10.0, 1.0), (100.0, 1.0)]
    kappa_pi_levels = [0.1, 0.3, 0.6]
    c2_levels = [0.2, 0.5, 1.0]
    tau_pi_levels = [5, 10, 40]
    reveal_ratio_levels = [0.15, 0.25, 0.45]
    warmup_iters_levels = [100, 500, 1000]
    unpin_levels = [true, false]
    K_init_levels = [3, 5, 10]
    wf_levels = ["db2", "db4", "db6", "la8"]
    
    # Create all combinations
    parameter_combinations = []
    exp_id = 1
    
    for alpha_prior in alpha_prior_levels
        for kappa_pi in kappa_pi_levels
            for c2 in c2_levels
                for tau_pi in tau_pi_levels
                    for reveal_ratio in reveal_ratio_levels
                        for warmup_iters in warmup_iters_levels
                            for unpin in unpin_levels
                                for K_init in K_init_levels
                                    for wf in wf_levels
                                        params = Dict(
                                            :alpha_prior => alpha_prior,
                                            :kappa_pi => kappa_pi,
                                            :c2 => c2,
                                            :tau_pi => tau_pi,
                                            :reveal_ratio => reveal_ratio,
                                            :warmup_iters => warmup_iters,
                                            :unpin => unpin,
                                            :K_init => K_init,
                                            :wf => wf
                                        )
                                        push!(parameter_combinations, (exp_id, params))
                                        exp_id += 1
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    return parameter_combinations
end

# Function to analyze results
function analyze_results(results_file::String)
    println("\nAnalyzing results from $results_file...")
    
    # Read all results
    results = TOML.parsefile(results_file)
    
    # Extract F1 scores and parameters
    f1_scores = Float64[]
    experiment_data = []
    
    for (key, data) in results
        if startswith(key, "experiment_")
            f1 = get(data, "f1_score", 0.0)
            push!(f1_scores, f1)
            push!(experiment_data, data)
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
    
    println("\n" * "="^80)
    println("FACTORIAL EXPERIMENT RESULTS")
    println("="^80)
    println("Total experiments: $(length(f1_scores))")
    println("Best F1 Score: $(round(best_f1, digits=4))")
    println("Mean F1 Score: $(round(mean(f1_scores), digits=4))")
    println("Std F1 Score: $(round(std(f1_scores), digits=4))")
    
    println("\nBest Parameterization:")
    println("-"^40)
    for (param, value) in best_experiment
        if param != "experiment_id" && param != "f1_score" && param != "precision" && 
           param != "recall" && param != "accuracy" && param != "ari" && 
           param != "n_clusters" && param != "normal_cluster" && 
           param != "tp" && param != "fp" && param != "fn" && param != "tn" && param != "error"
            println("  $param: $value")
        end
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
                "wf=$(exp_data["wf"]), " *
                "K_init=$(exp_data["K_init"])")
    end
    
    return best_experiment
end

# Main function
function main()
    # Setup
    data_dir = joinpath(@__DIR__, "..", "..", "..", "data", "ArrowHead")
    results_file = joinpath(@__DIR__, "factorial_results.toml")
    
    # Clear results file
    if isfile(results_file)
        rm(results_file)
    end
    
    # Create parameter grid
    parameter_combinations = create_parameter_grid()
    total_experiments = length(parameter_combinations)
    
    println("Factorial Experiment Setup")
    println("="^50)
    println("Total experiments: $total_experiments")
    println("Results file: $results_file")
    println("Data directory: $data_dir")
    
    # Check if data directory exists
    if !isdir(data_dir)
        println("ERROR: Data directory not found: $data_dir")
        return
    end
    
    # Run experiments
    println("\nStarting experiments...")
    start_time = time()
    
    # Use Threads for parallel execution
    Threads.@threads for (exp_id, params) in parameter_combinations
        run_single_experiment(exp_id, params, data_dir, results_file)
    end
    
    end_time = time()
    elapsed_time = end_time - start_time
    
    println("\nAll experiments completed!")
    println("Total time: $(round(elapsed_time, digits=2)) seconds")
    println("Average time per experiment: $(round(elapsed_time / total_experiments, digits=2)) seconds")
    
    # Analyze results
    analyze_results(results_file)
end

# Run if called directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
