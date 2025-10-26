module Plotting

using Plots
using Statistics: mean
using StatsBase: countmap

export plot_dataset_before_clustering, plot_dataset_after_clustering, plot_clustering_comparison

"""
    plot_dataset_before_clustering(series, labels, t; title="Dataset Before Clustering", save_path=nothing)

Plot the original dataset before clustering, showing all time series colored by their true labels.

# Arguments
- `series`: Vector of time series matrices
- `labels`: Vector of binary labels (0 for normal, 1 for anomaly)
- `t`: Time vector
- `title`: Plot title
- `save_path`: Optional path to save the plot

# Returns
- Plot object
"""
function plot_dataset_before_clustering(series, labels, t; title="Dataset Before Clustering", save_path=nothing)
    # Create subplots for each dimension
    n_dims = size(series[1], 2)
    plots = []
    
    for dim in 1:n_dims
        p = plot(title="Dimension $dim", xlabel="Time", ylabel="Value", 
                legend=:topright, size=(800, 400))
        
        # Plot normal samples (label = 0)
        normal_indices = findall(==(0), labels)
        for idx in normal_indices
            plot!(p, t, series[idx][:, dim], 
                  color=:blue, alpha=0.6, linewidth=1, 
                  label=idx == normal_indices[1] ? "Normal" : "")
        end
        
        # Plot anomaly samples (label = 1)
        anomaly_indices = findall(==(1), labels)
        for idx in anomaly_indices
            plot!(p, t, series[idx][:, dim], 
                  color=:red, alpha=0.8, linewidth=2, 
                  label=idx == anomaly_indices[1] ? "Anomaly" : "")
        end
        
        push!(plots, p)
    end
    
    # Combine plots
    if n_dims == 1
        final_plot = plots[1]
    else
        final_plot = plot(plots..., layout=(n_dims, 1), size=(800, 400*n_dims))
    end
    
    plot!(final_plot, title=title, plot_titlefontsize=14)
    
    if save_path !== nothing
        savefig(final_plot, save_path)
        println("Plot saved to: $save_path")
    end
    
    return final_plot
end

"""
    plot_dataset_after_clustering(series, clusters, t; title="Dataset After Clustering", save_path=nothing)

Plot the dataset after clustering, showing time series colored by their cluster assignments.

# Arguments
- `series`: Vector of time series matrices
- `clusters`: Vector of cluster assignments
- `t`: Time vector
- `title`: Plot title
- `save_path`: Optional path to save the plot

# Returns
- Plot object
"""
function plot_dataset_after_clustering(series, clusters, t; title="Dataset After Clustering", save_path=nothing)
    # Create subplots for each dimension
    n_dims = size(series[1], 2)
    plots = []
    
    # Get unique clusters and assign colors
    unique_clusters = sort(unique(clusters))
    colors = [:blue, :red, :green, :orange, :purple, :brown, :pink, :gray, :olive, :cyan]
    
    for dim in 1:n_dims
        p = plot(title="Dimension $dim", xlabel="Time", ylabel="Value", 
                legend=:topright, size=(800, 400))
        
        # Plot each cluster
        for (i, cluster) in enumerate(unique_clusters)
            cluster_indices = findall(==(cluster), clusters)
            color = colors[mod(i-1, length(colors)) + 1]
            
            for (j, idx) in enumerate(cluster_indices)
                plot!(p, t, series[idx][:, dim], 
                      color=color, alpha=0.6, linewidth=1,
                      label=j == 1 ? "Cluster $cluster" : "")
            end
        end
        
        push!(plots, p)
    end
    
    # Combine plots
    if n_dims == 1
        final_plot = plots[1]
    else
        final_plot = plot(plots..., layout=(n_dims, 1), size=(800, 400*n_dims))
    end
    
    plot!(final_plot, title=title, plot_titlefontsize=14)
    
    if save_path !== nothing
        savefig(final_plot, save_path)
        println("Plot saved to: $save_path")
    end
    
    return final_plot
end

"""
    plot_clustering_comparison(series, true_labels, clusters, t; save_path=nothing)

Create a side-by-side comparison of the dataset before and after clustering.

# Arguments
- `series`: Vector of time series matrices
- `true_labels`: Vector of true binary labels
- `clusters`: Vector of cluster assignments
- `t`: Time vector
- `save_path`: Optional path to save the plot

# Returns
- Plot object
"""
function plot_clustering_comparison(series, true_labels, clusters, t; save_path=nothing)
    n_dims = size(series[1], 2)
    
    # Create before and after plots for each dimension
    all_plots = []
    
    for dim in 1:n_dims
        # Before clustering plot
        p_before = plot(title="Before Clustering (Dim $dim)", xlabel="Time", ylabel="Value", 
                       legend=:topright, size=(400, 300))
        
        # Plot normal samples
        normal_indices = findall(==(0), true_labels)
        for idx in normal_indices
            plot!(p_before, t, series[idx][:, dim], 
                  color=:blue, alpha=0.6, linewidth=1, 
                  label=idx == normal_indices[1] ? "Normal" : "")
        end
        
        # Plot anomaly samples
        anomaly_indices = findall(==(1), true_labels)
        for idx in anomaly_indices
            plot!(p_before, t, series[idx][:, dim], 
                  color=:red, alpha=0.8, linewidth=2, 
                  label=idx == anomaly_indices[1] ? "Anomaly" : "")
        end
        
        # After clustering plot
        p_after = plot(title="After Clustering (Dim $dim)", xlabel="Time", ylabel="Value", 
                      legend=:topright, size=(400, 300))
        
        # Get unique clusters and assign colors
        unique_clusters = sort(unique(clusters))
        colors = [:blue, :red, :green, :orange, :purple, :brown, :pink, :gray, :olive, :cyan]
        
        for (i, cluster) in enumerate(unique_clusters)
            cluster_indices = findall(==(cluster), clusters)
            color = colors[mod(i-1, length(colors)) + 1]
            
            for (j, idx) in enumerate(cluster_indices)
                plot!(p_after, t, series[idx][:, dim], 
                      color=color, alpha=0.6, linewidth=1,
                      label=j == 1 ? "Cluster $cluster" : "")
            end
        end
        
        push!(all_plots, p_before, p_after)
    end
    
    # Create layout: 2 columns (before/after) × n_dims rows
    layout = @layout([grid(n_dims, 2)])
    final_plot = plot(all_plots..., layout=layout, size=(800, 300*n_dims))
    
    if save_path !== nothing
        savefig(final_plot, save_path)
        println("Comparison plot saved to: $save_path")
    end
    
    return final_plot
end

"""
    plot_clustering_summary(series, true_labels, clusters, t; save_dir=nothing)

Create a comprehensive summary of clustering results including before/after plots and statistics.

# Arguments
- `series`: Vector of time series matrices
- `true_labels`: Vector of true binary labels
- `clusters`: Vector of cluster assignments
- `t`: Time vector
- `save_dir`: Optional directory to save plots

# Returns
- Dictionary of plot objects
"""
function plot_clustering_summary(series, true_labels, clusters, t; save_dir=nothing)
    plots_dict = Dict()
    
    # Create plots directory if saving
    if save_dir !== nothing
        mkpath(save_dir)
    end
    
    # Plot 1: Before clustering
    before_path = save_dir !== nothing ? joinpath(save_dir, "dataset_before_clustering.png") : nothing
    plots_dict[:before] = plot_dataset_before_clustering(series, true_labels, t; 
                                                        title="Dataset Before Clustering",
                                                        save_path=before_path)
    
    # Plot 2: After clustering
    after_path = save_dir !== nothing ? joinpath(save_dir, "dataset_after_clustering.png") : nothing
    plots_dict[:after] = plot_dataset_after_clustering(series, clusters, t; 
                                                      title="Dataset After Clustering",
                                                      save_path=after_path)
    
    # Plot 3: Comparison
    comparison_path = save_dir !== nothing ? joinpath(save_dir, "clustering_comparison.png") : nothing
    plots_dict[:comparison] = plot_clustering_comparison(series, true_labels, clusters, t; 
                                                        save_path=comparison_path)
    
    # Print summary statistics
    println("\nClustering Summary:")
    println("==================")
    println("Total samples: $(length(series))")
    println("True anomalies: $(sum(true_labels))")
    println("True normal: $(sum(1 .- true_labels))")
    println("Number of clusters found: $(length(unique(clusters)))")
    println("Cluster sizes: $(countmap(clusters))")
    
    return plots_dict
end

end # module
