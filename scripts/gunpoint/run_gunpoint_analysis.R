#!/usr/bin/env Rscript

# GunPoint Dataset Analysis with WICMAD
# This script loads the GunPoint dataset, runs WICMAD clustering, and visualizes results

# Set working directory to scripts folder
# setwd("scripts")

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(cluster)
library(mclust)  # For adjustedRandIndex
library(numDeriv)

# Load WICMAD package
library(devtools)
devtools::load_all()

# Function to load GunPoint dataset with imbalanced sampling
load_gunpoint_data <- function(data_dir) {
  cat("Loading GunPoint dataset...\n")
  
  # Load training and test data
  train_file <- file.path(data_dir, "GunPoint_TRAIN.txt")
  test_file <- file.path(data_dir, "GunPoint_TEST.txt")
  
  if (!file.exists(train_file) || !file.exists(test_file)) {
    stop("GunPoint data files not found in:", data_dir)
  }
  
  # Read training data
  train_data <- read.table(train_file, header = FALSE)
  test_data <- read.table(test_file, header = FALSE)
  
  # Combine train and test data
  combined_data <- rbind(train_data, test_data)
  
  cat("Total samples:", nrow(combined_data), "\n")
  cat("Training samples:", nrow(train_data), "\n")
  cat("Test samples:", nrow(test_data), "\n")
  
  # Use all observations in the normal class, then make final dataset 5% anomalies
  # (This is already handled in the imbalanced sampling below)
  
  # Extract labels (first column) and time series (remaining columns)
  labels <- combined_data[, 1]
  time_series <- as.matrix(combined_data[, -1])
  
  # Convert labels: 1.0 = Gun-Draw (normal), 2.0 = Point (anomaly)
  # We'll treat Gun-Draw as normal (0) and Point as anomaly (1)
  binary_labels <- ifelse(labels == 1.0, 0, 1)
  
  cat("Original class distribution:\n")
  print(table(labels))
  cat("Binary labels (0=Gun-Draw/Normal, 1=Point/Anomaly):\n")
  print(table(binary_labels))
  
  # Create imbalanced dataset with 5% anomalies
  normal_indices <- which(binary_labels == 0)
  anomaly_indices <- which(binary_labels == 1)
  
  # Calculate how many anomalies we need for 5% of total
  total_samples <- length(binary_labels)
  n_anomalies_needed <- max(1, round(total_samples * 0.05))
  n_normal_needed <- total_samples - n_anomalies_needed
  
  # Sample the required number of normal and anomaly samples
  if (length(normal_indices) >= n_normal_needed) {
    selected_normal <- sample(normal_indices, n_normal_needed)
  } else {
    selected_normal <- normal_indices
  }
  
  if (length(anomaly_indices) >= n_anomalies_needed) {
    selected_anomaly <- sample(anomaly_indices, n_anomalies_needed)
  } else {
    selected_anomaly <- anomaly_indices
  }
  
  # Combine selected indices
  selected_indices <- c(selected_normal, selected_anomaly)
  
  # Create imbalanced dataset
  imbalanced_time_series <- time_series[selected_indices, ]
  imbalanced_labels <- binary_labels[selected_indices]
  
  # Shuffle the data
  shuffle_indices <- sample(length(imbalanced_labels))
  imbalanced_time_series <- imbalanced_time_series[shuffle_indices, ]
  imbalanced_labels <- imbalanced_labels[shuffle_indices]
  
  cat("Final imbalanced dataset:\n")
  cat("Total samples:", nrow(imbalanced_time_series), "\n")
  cat("Normal (Gun-Draw):", sum(imbalanced_labels == 0), "\n")
  cat("Anomaly (Point):", sum(imbalanced_labels == 1), "\n")
  cat("Anomaly percentage:", round(mean(imbalanced_labels == 1) * 100, 1), "%\n")
  
  # Set target dimension to 16 for testing
  target_dim <- 16
  cat("Using target dimension for testing:", target_dim, "\n")
  
  return(list(
    train_data = imbalanced_time_series,
    train_labels = imbalanced_labels,
    majority_class = 0,  # Gun-Draw is the majority class
    target_dim = target_dim
  ))
}

# Function to interpolate time series to target dimension
interpolate_series <- function(series, target_dim = 16) {
  # Handle edge cases first
  if (length(series) < 2) {
    return(rep(series[1], target_dim))
  }
  
  # Check for NA values
  if (any(is.na(series))) {
    series <- na.omit(series)
    if (length(series) < 2) {
      return(rep(series[1], target_dim))
    }
  }
  
  # If already correct length, return as is
  if (length(series) == target_dim) {
    return(series)
  }
  
  # Use spline interpolation which is more robust
  tryCatch({
    x_old <- seq(0, 1, length.out = length(series))
    x_new <- seq(0, 1, length.out = target_dim)
    interpolated <- spline(x_old, series, xout = x_new)$y
    return(interpolated)
  }, error = function(e) {
    # Fallback: simple resampling
    indices <- round(seq(1, length(series), length.out = target_dim))
    return(series[indices])
  })
}

# Function to prepare data for WICMAD (raw time series)
prepare_wicmad_data <- function(series_data, labels, target_dim = NULL) {
  cat("Preparing WICMAD data for raw time series...\n")
  
  # Find the highest power of 2 supported by the data
  if (is.null(target_dim)) {
    if (is.matrix(series_data)) {
      max_dim <- ncol(series_data)
    } else {
      max_dim <- max(sapply(series_data, length))
    }
    target_dim <- 2^floor(log2(max_dim))
    cat("Using p parameter (highest power of 2):", target_dim, "\n")
  }
  
  # Handle both matrix and list inputs
  if (is.matrix(series_data)) {
    # Convert matrix to list of series (each row is a series)
    series_list <- lapply(seq_len(nrow(series_data)), function(i) {
      as.numeric(series_data[i, ])
    })
  } else {
    # Already a list
    series_list <- series_data
  }
  
  # Interpolate each series to target dimension
  interpolated_series <- lapply(series_list, function(series) {
    interpolate_series(series, target_dim)
  })
  
  # Convert to matrix format for WICMAD
  Y_matrix <- do.call(rbind, interpolated_series)
  
  # Create time grid
  t <- seq(0, 1, length.out = target_dim)
  
  # Convert to list format for WICMAD (each row becomes a P x 1 matrix)
  Y_list <- lapply(seq_len(nrow(Y_matrix)), function(i) {
    matrix(Y_matrix[i, ], nrow = target_dim, ncol = 1)
  })
  
  cat("WICMAD data prepared:\n")
  cat("- Number of series:", length(Y_list), "\n")
  cat("- Time points per series:", target_dim, "\n")
  cat("- Dimensions: P=", target_dim, "(time points), M=1 (univariate)\n")
  
  return(list(
    Y = Y_list,
    t = t,
    labels = labels
  ))
}

# Function to prepare derivative data for WICMAD
prepare_derivative_data <- function(series_data, labels, target_dim = NULL) {
  cat("Preparing derivative data for WICMAD...\n")
  
  # Find the highest power of 2 supported by the data
  if (is.null(target_dim)) {
    if (is.matrix(series_data)) {
      max_dim <- ncol(series_data)
    } else {
      max_dim <- max(sapply(series_data, length))
    }
    target_dim <- 2^floor(log2(max_dim))
    cat("Using p parameter for derivatives (highest power of 2):", target_dim, "\n")
  }
  
  # Handle both matrix and list inputs
  if (is.matrix(series_data)) {
    # Convert matrix to list of series (each row is a series)
    series_list <- lapply(seq_len(nrow(series_data)), function(i) {
      as.numeric(series_data[i, ])
    })
  } else {
    # Already a list
    series_list <- series_data
  }
  
  # Interpolate each series to target dimension
  interpolated_series <- lapply(series_list, function(series) {
    interpolate_series(series, target_dim)
  })
  
  # Convert to matrix format
  Y_matrix <- do.call(rbind, interpolated_series)
  
  # Create time grid
  t <- seq(0, 1, length.out = target_dim)
  
  # Calculate first and second derivatives using numDeriv
  first_deriv <- t(apply(Y_matrix, 1, function(row) {
    # Create interpolation function for smooth derivatives
    interp_func <- splinefun(t, row, method = "natural")
    
    # Compute first derivative using numDeriv
    sapply(t, function(x) {
      tryCatch({
        grad(interp_func, x)
      }, error = function(e) {
        # Fallback to finite differences if grad fails
        if (x == 0) return((interp_func(t[2]) - interp_func(t[1])) / (t[2] - t[1]))
        if (x == 1) return((interp_func(t[target_dim]) - interp_func(t[target_dim-1])) / (t[target_dim] - t[target_dim-1]))
        h <- (t[2] - t[1])
        return((interp_func(x + h) - interp_func(x - h)) / (2 * h))
      })
    })
  }))
  
  second_deriv <- t(apply(Y_matrix, 1, function(row) {
    # Create interpolation function for smooth derivatives
    interp_func <- splinefun(t, row, method = "natural")
    
    # Compute second derivative using numDeriv
    sapply(t, function(x) {
      tryCatch({
        hessian(interp_func, x)
      }, error = function(e) {
        # Fallback to finite differences if hessian fails
        if (x == 0 || x == 1) return(0)
        h <- (t[2] - t[1])
        return((interp_func(x + h) - 2*interp_func(x) + interp_func(x - h)) / (h^2))
      })
    })
  }))
  
  # Combine original, first derivative, and second derivative
  combined_data <- cbind(Y_matrix, first_deriv, second_deriv)
  
  # Convert to list format for WICMAD (each row becomes a P x 3 matrix)
  Y_list <- lapply(seq_len(nrow(combined_data)), function(i) {
    matrix(combined_data[i, ], nrow = target_dim, ncol = 3)
  })
  
  cat("Derivative data prepared:\n")
  cat("- Number of series:", length(Y_list), "\n")
  cat("- Time points per series:", target_dim, "\n")
  cat("- Dimensions: P=", target_dim, "(time points), M=3 (original + 1st + 2nd derivatives)\n")
  
  return(list(
    Y = Y_list,
    t = t,
    labels = labels,
    original_data = Y_matrix,
    first_deriv_data = first_deriv,
    second_deriv_data = second_deriv
  ))
}


# Function to run WICMAD analysis
run_wicmad_analysis <- function(wicmad_data, analysis_name) {
  cat("\n=== Running WICMAD analysis:", analysis_name, "===\n")
  
  # Run WICMAD
  wicmad_result <- wicmad(
    Y = wicmad_data$Y,
    t = wicmad_data$t,
    n_iter = 8000,
    burn = 3000,
    warmup_iters = 500,
    unpin = FALSE
  )
  
  # Extract cluster assignments
  cluster_assignments <- wicmad_result$Z[nrow(wicmad_result$Z), ]
  
  cat("WICMAD completed. Number of clusters found:", length(unique(cluster_assignments)), "\n")
  
  # Map WICMAD clusters to binary classification (normal vs anomaly)
  normal_samples <- which(wicmad_data$labels == 0)
  cluster_counts <- table(cluster_assignments[normal_samples])
  normal_cluster <- as.numeric(names(cluster_counts)[which.max(cluster_counts)])
  binary_cluster_assignments <- ifelse(cluster_assignments == normal_cluster, 0, 1)
  
  # Calculate metrics
  true_labels <- wicmad_data$labels
  pred_labels <- binary_cluster_assignments
  
  # Confusion matrix
  cm <- table(true_labels, pred_labels)
  cat("Confusion Matrix:\n")
  print(cm)
  
  # Calculate metrics
  if (nrow(cm) == 2 && ncol(cm) == 2) {
    tp <- cm[2, 2]  # True positives (anomaly correctly identified)
    tn <- cm[1, 1]  # True negatives (normal correctly identified)
    fp <- cm[1, 2]  # False positives (normal misclassified as anomaly)
    fn <- cm[2, 1]  # False negatives (anomaly misclassified as normal)
    
    accuracy <- (tp + tn) / (tp + tn + fp + fn)
    precision <- tp / (tp + fp)
    recall <- tp / (tp + fn)
    f1 <- 2 * (precision * recall) / (precision + recall)
    
    cat("Metrics:\n")
    cat("Accuracy:", round(accuracy, 4), "\n")
    cat("Precision:", round(precision, 4), "\n")
    cat("Recall:", round(recall, 4), "\n")
    cat("F1-Score:", round(f1, 4), "\n")
  }
  
  # Adjusted Rand Index
  ari <- adjustedRandIndex(true_labels, pred_labels)
  cat("Adjusted Rand Index:", round(ari, 4), "\n")
  
  return(list(
    cluster_assignments = cluster_assignments,
    binary_cluster_assignments = binary_cluster_assignments,
    confusion_matrix = cm,
    accuracy = if (exists("accuracy")) accuracy else NA,
    precision = if (exists("precision")) precision else NA,
    recall = if (exists("recall")) recall else NA,
    f1 = if (exists("f1")) f1 else NA,
    ari = ari
  ))
}

# Function to plot overlapped data
plot_overlapped_data <- function(data, labels, title = "GunPoint Dataset", max_series = NULL) {
  cat("Creating overlapped data plot...\n")
  
  # Use all series if max_series is NULL
  if (is.null(max_series)) {
    series_indices <- seq_len(nrow(data))
  } else {
    n_series <- min(max_series, nrow(data))
    series_indices <- sample(seq_len(nrow(data)), n_series)
  }
  
  # Prepare data for plotting
  plot_data <- data.frame()
  
  for (i in series_indices) {
    series_df <- data.frame(
      time = seq_len(ncol(data)),
      value = as.numeric(data[i, ]),
      series_id = paste("Series", i),
      label = labels[i],
      label_name = ifelse(labels[i] == 0, "Gun-Draw (Normal)", "Point (Anomaly)")
    )
    plot_data <- rbind(plot_data, series_df)
  }
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = .data$time, y = .data$value, color = .data$label_name, group = .data$series_id)) +
    geom_line(alpha = 0.7, size = 0.5) +
    labs(title = title, x = "Time", y = "Value", color = "Class") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(p)
}

# Function to plot clustering results with cluster assignments
plot_clustering_results <- function(data, true_labels, cluster_assignments, title = "GunPoint Dataset", max_series = NULL) {
  cat("Creating clustering results plot...\n")
  
  # Use all series if max_series is NULL
  if (is.null(max_series)) {
    series_indices <- seq_len(nrow(data))
  } else {
    n_series <- min(max_series, nrow(data))
    series_indices <- sample(seq_len(nrow(data)), n_series)
  }
  
  # Prepare data for plotting
  plot_data <- data.frame()
  
  for (i in series_indices) {
    series_df <- data.frame(
      time = seq_len(ncol(data)),
      value = as.numeric(data[i, ]),
      series_id = paste("Series", i),
      true_label = true_labels[i],
      cluster = cluster_assignments[i],
      cluster_name = ifelse(cluster_assignments[i] == 0, "Normal", "Anomaly")
    )
    plot_data <- rbind(plot_data, series_df)
  }
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = .data$time, y = .data$value, color = .data$cluster_name, group = .data$series_id)) +
    geom_line(alpha = 0.7, size = 0.5) +
    labs(title = title, x = "Time", y = "Value", color = "Class") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(p)
}

# Function to plot derivative clustering results (3 separate channels)
plot_derivative_clustering_results <- function(original_data, first_deriv_data, second_deriv_data, true_labels, cluster_assignments, title = "GunPoint - After Clustering", max_series = NULL) {
  cat("Creating derivative clustering results plot...\n")
  
  # Use all series if max_series is NULL
  if (is.null(max_series)) {
    series_indices <- seq_len(nrow(original_data))
  } else {
    n_series <- min(max_series, nrow(original_data))
    series_indices <- sample(seq_len(nrow(original_data)), n_series)
  }
  
  # Prepare data for plotting
  plot_data <- data.frame()
  
  for (i in series_indices) {
    # Original data
    orig_df <- data.frame(
      time = seq_len(ncol(original_data)),
      value = as.numeric(original_data[i, ]),
      series_id = paste("Series", i),
      channel = "Original",
      cluster = cluster_assignments[i],
      cluster_name = ifelse(cluster_assignments[i] == 0, "Normal", "Anomaly")
    )
    
    # First derivative
    first_df <- data.frame(
      time = seq_len(ncol(first_deriv_data)),
      value = as.numeric(first_deriv_data[i, ]),
      series_id = paste("Series", i),
      channel = "First Derivative",
      cluster = cluster_assignments[i],
      cluster_name = ifelse(cluster_assignments[i] == 0, "Normal", "Anomaly")
    )
    
    # Second derivative
    second_df <- data.frame(
      time = seq_len(ncol(second_deriv_data)),
      value = as.numeric(second_deriv_data[i, ]),
      series_id = paste("Series", i),
      channel = "Second Derivative",
      cluster = cluster_assignments[i],
      cluster_name = ifelse(cluster_assignments[i] == 0, "Normal", "Anomaly")
    )
    
    plot_data <- rbind(plot_data, orig_df, first_df, second_df)
  }
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = .data$time, y = .data$value, color = .data$cluster_name, group = .data$series_id)) +
    geom_line(alpha = 0.7, size = 0.5) +
    facet_wrap(~ channel, scales = "free_y", ncol = 1) +
    labs(title = title, x = "Time", y = "Value", color = "Class") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(p)
}


# Main analysis function
main <- function() {
  cat("=== GunPoint Dataset Analysis with WICMAD ===\n")
  
  # Load data
  cat("\n1. Loading GunPoint dataset...\n")
  data <- load_gunpoint_data("../data/GunPoint")
  
  cat("Label distribution:\n")
  print(table(data$train_labels))
  cat("Majority class used as normal:", data$majority_class, "\n")
  cat("First few labels:", head(data$train_labels, 10), "\n")
  
  # Plot original data (overlapped)
  cat("\n2. Creating original data visualization...\n")
  original_plot <- plot_overlapped_data(data$train_data, data$train_labels, 
                                       "GunPoint Dataset - Before Clustering")
  
  # Save original data plot
  pdf("../plots/gunpoint/gunpoint_imbalanced_original_data.pdf", width = 12, height = 8)
  print(original_plot)
  dev.off()
  cat("Imbalanced original data plot saved to ../plots/gunpoint/gunpoint_imbalanced_original_data.pdf\n")
  
  # Run analysis on raw time series
  cat("\n3. Analyzing raw time series data...\n")
  raw_wicmad_data <- prepare_wicmad_data(data$train_data, data$train_labels, data$target_dim)
  raw_results <- run_wicmad_analysis(raw_wicmad_data, "Raw Time Series")
  
  # Plot raw time series clustering results
  raw_clustering_plot <- plot_clustering_results(
    data$train_data, 
    data$train_labels, 
    raw_results$cluster_assignments,
    "GunPoint Dataset - After Clustering"
  )
  
  pdf("../plots/gunpoint/gunpoint_raw_curves_clustering.pdf", width = 12, height = 8)
  print(raw_clustering_plot)
  dev.off()
  cat("Raw curves clustering plot saved to ../plots/gunpoint/gunpoint_raw_curves_clustering.pdf\n")
  
  # Run analysis on derivatives
  cat("\n4. Analyzing derivative data...\n")
  deriv_wicmad_data <- prepare_derivative_data(data$train_data, data$train_labels, data$target_dim)
  deriv_results <- run_wicmad_analysis(deriv_wicmad_data, "Derivatives")
  
  # Plot derivative clustering results (3 separate channels)
  deriv_data_prep <- prepare_derivative_data(data$train_data, data$train_labels, data$target_dim)
  deriv_clustering_plot <- plot_derivative_clustering_results(
    deriv_data_prep$original_data,
    deriv_data_prep$first_deriv_data,
    deriv_data_prep$second_deriv_data,
    deriv_results$cluster_assignments, 
    deriv_results$cluster_assignments,
    "GunPoint Dataset - After Clustering"
  )
  
  pdf("../plots/gunpoint/gunpoint_derivatives_clustering.pdf", width = 12, height = 8)
  print(deriv_clustering_plot)
  dev.off()
  cat("Derivatives clustering plot saved to ../plots/gunpoint/gunpoint_derivatives_clustering.pdf\n")
  
  
  # Summary
  cat("\n=== Analysis Summary ===\n")
  cat("Raw Time Series - ARI:", round(raw_results$ari, 4), "\n")
  cat("Derivatives - ARI:", round(deriv_results$ari, 4), "\n")
  
  cat("\nAnalysis completed successfully!\n")
}

main()
