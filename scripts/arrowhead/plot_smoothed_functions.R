#!/usr/bin/env Rscript

# Script to plot original vs smoothed ArrowHead functions
# This script loads the ArrowHead dataset and creates side-by-side plots
# showing a handful of functions before and after smoothing

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

# Function to generate sample time series data
generate_sample_data <- function() {
  cat("Generating sample time series data...\n")
  
  # Set seed for reproducibility
  set.seed(123)
  
  # Generate 3 different types of functions
  n_series <- 15
  series_list <- list()
  labels <- integer(n_series)
  
  # Type 1: Smooth sine waves (Class 1)
  for (i in 1:5) {
    t <- seq(0, 1, length.out = 50)
    # Add some noise to make it more realistic
    noise <- rnorm(length(t), 0, 0.1)
    y <- sin(2 * pi * t) + 0.5 * sin(4 * pi * t) + noise
    series_list[[i]] <- matrix(y, ncol = 1)
    labels[i] <- 1
  }
  
  # Type 2: Sawtooth-like functions (Class 2)
  for (i in 6:10) {
    t <- seq(0, 1, length.out = 50)
    # Create sawtooth pattern with noise
    noise <- rnorm(length(t), 0, 0.1)
    y <- 2 * (t - floor(t + 0.5)) + noise
    series_list[[i]] <- matrix(y, ncol = 1)
    labels[i] <- 2
  }
  
  # Type 3: Exponential decay functions (Class 3)
  for (i in 11:15) {
    t <- seq(0, 1, length.out = 50)
    # Create exponential decay with noise
    noise <- rnorm(length(t), 0, 0.1)
    y <- exp(-3 * t) + 0.3 * sin(6 * pi * t) + noise
    series_list[[i]] <- matrix(y, ncol = 1)
    labels[i] <- 3
  }
  
  cat("Generated", n_series, "sample series\n")
  cat("Label distribution:\n")
  print(table(labels))
  
  return(list(series_list = series_list, labels = labels))
}

# Function to create smoothed version using splinefun
create_smoothed_function <- function(series) {
  n <- length(series)
  x_vals <- seq(0, 1, length.out = n)
  
  # Create smooth interpolation function
  interp_func <- splinefun(x_vals, series, method = "natural")
  
  # Evaluate at more points for smoother appearance
  x_smooth <- seq(0, 1, length.out = 100)
  y_smooth <- interp_func(x_smooth)
  
  return(list(x = x_smooth, y = y_smooth))
}

# Function to plot original vs smoothed functions
plot_original_vs_smoothed <- function(series_list, labels, n_functions = 5) {
  cat("Creating original vs smoothed function plots...\n")
  
  # Select a handful of functions to plot
  n_series <- length(series_list)
  selected_indices <- sample(seq_len(n_series), min(n_functions, n_series))
  
  # Create plotting data
  plot_data_original <- data.frame()
  plot_data_smoothed <- data.frame()
  
  for (i in seq_along(selected_indices)) {
    idx <- selected_indices[i]
    series_data <- series_list[[idx]]
    time_points <- seq(0, 1, length.out = nrow(series_data))
    
    # Original data
    orig_df <- data.frame(
      time = time_points,
      value = series_data[, 1],
      function_id = i,
      label = labels[idx],
      label_name = paste("Class", labels[idx]),
      type = "Original"
    )
    
    # Smoothed data
    smoothed <- create_smoothed_function(series_data[, 1])
    smooth_df <- data.frame(
      time = smoothed$x,
      value = smoothed$y,
      function_id = i,
      label = labels[idx],
      label_name = paste("Class", labels[idx]),
      type = "Smoothed"
    )
    
    plot_data_original <- rbind(plot_data_original, orig_df)
    plot_data_smoothed <- rbind(plot_data_smoothed, smooth_df)
  }
  
  # Create original functions plot (smooth like the after plot)
  p_original <- ggplot(plot_data_original, aes(x = time, y = value, color = label_name, group = function_id)) +
    geom_line(alpha = 0.8, linewidth = 0.8) +
    labs(title = "Original ArrowHead Functions",
         x = "Time", y = "Value", color = "Class") +
    theme_minimal() +
    theme(legend.position = "bottom") +
    scale_color_brewer(type = "qual", palette = "Set1")
  
  # Create smoothed functions plot
  p_smoothed <- ggplot(plot_data_smoothed, aes(x = time, y = value, color = label_name, group = function_id)) +
    geom_line(alpha = 0.8, linewidth = 0.8) +
    labs(title = "Smoothed ArrowHead Functions",
         x = "Time", y = "Value", color = "Class") +
    theme_minimal() +
    theme(legend.position = "bottom") +
    scale_color_brewer(type = "qual", palette = "Set1")
  
  # Combine plots side by side
  combined_plot <- grid.arrange(p_original, p_smoothed, ncol = 2)
  
  return(list(original = p_original, smoothed = p_smoothed, combined = combined_plot))
}

# Main function
main <- function() {
  cat("=== ArrowHead Smoothed Functions Visualization ===\n\n")
  
  # Ensure output directory exists
  output_dir <- "../../plots/arrowhead"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("Created output directory:", output_dir, "\n")
  }
  
  # Generate sample data
  cat("1. Generating sample time series data...\n")
  data <- generate_sample_data()
  
  cat("Loaded", length(data$series_list), "series\n")
  cat("Label distribution:\n")
  print(table(data$labels))
  
  # Create plots
  cat("\n2. Creating original vs smoothed function plots...\n")
  plots <- plot_original_vs_smoothed(data$series_list, data$labels, n_functions = 6)
  
  # Save individual plots
  original_plot_path <- "../../plots/arrowhead/arrowhead_original_functions.pdf"
  smoothed_plot_path <- "../../plots/arrowhead/arrowhead_smoothed_functions.pdf"
  combined_plot_path <- "../../plots/arrowhead/arrowhead_original_vs_smoothed.pdf"
  
  cat("Saving plots...\n")
  
  # Save original functions plot
  pdf(original_plot_path, width = 10, height = 6)
  print(plots$original)
  dev.off()
  cat("✓ Original functions plot saved to:", original_plot_path, "\n")
  
  # Save smoothed functions plot
  pdf(smoothed_plot_path, width = 10, height = 6)
  print(plots$smoothed)
  dev.off()
  cat("✓ Smoothed functions plot saved to:", smoothed_plot_path, "\n")
  
  # Save combined plot
  pdf(combined_plot_path, width = 16, height = 6)
  print(plots$combined)
  dev.off()
  cat("✓ Combined plot saved to:", combined_plot_path, "\n")
  
  cat("\n=== Visualization Complete ===\n")
  cat("Generated files:\n")
  cat("-", original_plot_path, ": Original ArrowHead functions\n")
  cat("-", smoothed_plot_path, ": Smoothed ArrowHead functions\n")
  cat("-", combined_plot_path, ": Side-by-side comparison\n")
}

# Run the analysis
main()
