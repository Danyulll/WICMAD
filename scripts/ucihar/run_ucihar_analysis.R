#!/usr/bin/env Rscript

# UCI HAR Dataset Analysis with WICMAD
# This script loads the UCI HAR dataset, runs WICMAD clustering, and visualizes results

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(cluster)
library(mclust)  # For adjustedRandIndex

# Load WICMAD package
library(devtools)
devtools::load_all()

# Function to load UCI HAR dataset and create imbalanced data
load_ucihar_data <- function(data_dir) {
  cat("Loading UCI HAR dataset...\n")
  
  # Load activity labels
  activity_labels <- read.table(file.path(data_dir, "UCI HAR Dataset", "UCI HAR Dataset", "activity_labels.txt"), 
                                col.names = c("id", "activity"))
  
  # Load training data
  cat("Loading training data...\n")
  train_labels <- read.table(file.path(data_dir, "UCI HAR Dataset", "UCI HAR Dataset", "train", "y_train.txt"), 
                             col.names = c("activity_id"))
  
  # Load test data
  cat("Loading test data...\n")
  test_labels <- read.table(file.path(data_dir, "UCI HAR Dataset", "UCI HAR Dataset", "test", "y_test.txt"), 
                            col.names = c("activity_id"))
  
  # Load inertial signals for training
  cat("Loading training inertial signals...\n")
  train_signals <- list(
    body_acc_x = read.table(file.path(data_dir, "UCI HAR Dataset", "UCI HAR Dataset", "train", "Inertial Signals", "body_acc_x_train.txt")),
    body_acc_y = read.table(file.path(data_dir, "UCI HAR Dataset", "UCI HAR Dataset", "train", "Inertial Signals", "body_acc_y_train.txt")),
    body_acc_z = read.table(file.path(data_dir, "UCI HAR Dataset", "UCI HAR Dataset", "train", "Inertial Signals", "body_acc_z_train.txt")),
    body_gyro_x = read.table(file.path(data_dir, "UCI HAR Dataset", "UCI HAR Dataset", "train", "Inertial Signals", "body_gyro_x_train.txt")),
    body_gyro_y = read.table(file.path(data_dir, "UCI HAR Dataset", "UCI HAR Dataset", "train", "Inertial Signals", "body_gyro_y_train.txt")),
    body_gyro_z = read.table(file.path(data_dir, "UCI HAR Dataset", "UCI HAR Dataset", "train", "Inertial Signals", "body_gyro_z_train.txt"))
  )
  
  # Load inertial signals for test
  cat("Loading test inertial signals...\n")
  test_signals <- list(
    body_acc_x = read.table(file.path(data_dir, "UCI HAR Dataset", "UCI HAR Dataset", "test", "Inertial Signals", "body_acc_x_test.txt")),
    body_acc_y = read.table(file.path(data_dir, "UCI HAR Dataset", "UCI HAR Dataset", "test", "Inertial Signals", "body_acc_y_test.txt")),
    body_acc_z = read.table(file.path(data_dir, "UCI HAR Dataset", "UCI HAR Dataset", "test", "Inertial Signals", "body_acc_z_test.txt")),
    body_gyro_x = read.table(file.path(data_dir, "UCI HAR Dataset", "UCI HAR Dataset", "test", "Inertial Signals", "body_gyro_x_test.txt")),
    body_gyro_y = read.table(file.path(data_dir, "UCI HAR Dataset", "UCI HAR Dataset", "test", "Inertial Signals", "body_gyro_y_test.txt")),
    body_gyro_z = read.table(file.path(data_dir, "UCI HAR Dataset", "UCI HAR Dataset", "test", "Inertial Signals", "body_gyro_z_test.txt"))
  )
  
  # Create multivariate time series for training data
  cat("Creating multivariate time series for training data...\n")
  train_series <- list()
  for (i in seq_len(nrow(train_labels))) {
    # Combine 6 sensor signals into one matrix (128 time points x 6 sensors)
    series_matrix <- matrix(0, nrow = 128, ncol = 6)
    series_matrix[, 1] <- as.numeric(train_signals$body_acc_x[i, ])
    series_matrix[, 2] <- as.numeric(train_signals$body_acc_y[i, ])
    series_matrix[, 3] <- as.numeric(train_signals$body_acc_z[i, ])
    series_matrix[, 4] <- as.numeric(train_signals$body_gyro_x[i, ])
    series_matrix[, 5] <- as.numeric(train_signals$body_gyro_y[i, ])
    series_matrix[, 6] <- as.numeric(train_signals$body_gyro_z[i, ])
    
    train_series[[i]] <- series_matrix
  }
  
  # Create multivariate time series for test data
  cat("Creating multivariate time series for test data...\n")
  test_series <- list()
  for (i in seq_len(nrow(test_labels))) {
    # Combine 6 sensor signals into one matrix (128 time points x 6 sensors)
    series_matrix <- matrix(0, nrow = 128, ncol = 6)
    series_matrix[, 1] <- as.numeric(test_signals$body_acc_x[i, ])
    series_matrix[, 2] <- as.numeric(test_signals$body_acc_y[i, ])
    series_matrix[, 3] <- as.numeric(test_signals$body_acc_z[i, ])
    series_matrix[, 4] <- as.numeric(test_signals$body_gyro_x[i, ])
    series_matrix[, 5] <- as.numeric(test_signals$body_gyro_y[i, ])
    series_matrix[, 6] <- as.numeric(test_signals$body_gyro_z[i, ])
    
    test_series[[i]] <- series_matrix
  }
  
  # Combine train and test data
  all_series <- c(train_series, test_series)
  all_labels <- c(train_labels$activity_id, test_labels$activity_id)
  
  cat("Dataset dimensions:\n")
  cat("Number of samples:", length(all_series), "\n")
  cat("Time points per sample:", nrow(all_series[[1]]), "\n")
  cat("Number of sensors:", ncol(all_series[[1]]), "\n")
  
  # Find the majority class to use as normal
  class_counts <- table(all_labels)
  majority_class <- as.numeric(names(class_counts)[which.max(class_counts)])
  minority_classes <- as.numeric(names(class_counts)[class_counts != max(class_counts)])
  
  cat("Class distribution:\n")
  print(class_counts)
  cat("Majority class (normal):", majority_class, "\n")
  cat("Minority classes (anomalies):", paste(minority_classes, collapse = ", "), "\n")
  
  # Convert to anomaly detection format (majority class = normal, others = anomaly)
  anomaly_labels <- ifelse(all_labels == majority_class, 0, 1)
  
  # Create imbalanced dataset with approximately 10% anomalies
  normal_indices <- which(anomaly_labels == 0)
  anomaly_indices <- which(anomaly_labels == 1)
  
  # Use all available normal samples
  selected_normal <- normal_indices
  n_normal_used <- length(selected_normal)
  
  # Calculate how many anomalies we need for 10% of the total dataset
  # Total dataset will be: n_normal_used + n_anomalies_needed
  # We want: n_anomalies_needed / (n_normal_used + n_anomalies_needed) = 0.10
  # Solving: n_anomalies_needed = 0.10 * (n_normal_used + n_anomalies_needed)
  # n_anomalies_needed = 0.10 * n_normal_used / (1 - 0.10) = 0.10 * n_normal_used / 0.90
  n_anomalies_needed <- max(1, round(0.10 * n_normal_used / 0.90))
  
  # Ensure we don't exceed available anomaly samples
  n_anomalies_needed <- min(n_anomalies_needed, length(anomaly_indices))
  
  # Sample the required number of anomaly samples
  selected_anomaly <- sample(anomaly_indices, n_anomalies_needed)
  
  # Combine selected indices
  selected_indices <- c(selected_normal, selected_anomaly)
  
  # Create imbalanced dataset
  imbalanced_series <- all_series[selected_indices]
  imbalanced_labels <- anomaly_labels[selected_indices]
  
  # Shuffle the data
  shuffle_indices <- sample(length(imbalanced_labels))
  imbalanced_series <- imbalanced_series[shuffle_indices]
  imbalanced_labels <- imbalanced_labels[shuffle_indices]
  
  cat("Final imbalanced dataset:\n")
  cat("Total samples:", length(imbalanced_series), "\n")
  cat("Normal (Class", majority_class, "):", sum(imbalanced_labels == 0), "\n")
  cat("Anomaly (Classes", paste(minority_classes, collapse = ", "), "):", sum(imbalanced_labels == 1), "\n")
  cat("Anomaly percentage:", round(mean(imbalanced_labels == 1) * 100, 1), "%\n")
  
  # Debug output for anomaly percentage
  cat("Debug: n_normal_used =", n_normal_used, ", n_anomalies_needed =", n_anomalies_needed, "\n")
  cat("Debug: Actual anomaly percentage =", round(mean(imbalanced_labels == 1) * 100, 1), "%\n")
  
  return(list(
    train_series = imbalanced_series,
    train_labels = imbalanced_labels,
    test_series = NULL,
    test_labels = NULL
  ))
}

# Function to interpolate time series to target dimension
interpolate_series <- function(series, target_dim = NULL) {
  # Find the highest power of 2 supported by the data
  if (is.null(target_dim)) {
    max_dim <- nrow(series)
    target_dim <- 2^floor(log2(max_dim))
    cat("Using p parameter (highest power of 2):", target_dim, "\n")
  }
  n <- nrow(series)
  if (n == target_dim) {
    return(series)
  }
  
  # Create time points for interpolation
  x_old <- seq(0, 1, length.out = n)
  x_new <- seq(0, 1, length.out = target_dim)
  
  # Interpolate each dimension
  interpolated <- matrix(0, nrow = target_dim, ncol = ncol(series))
  for (j in seq_len(ncol(series))) {
    interpolated[, j] <- approx(x_old, series[, j], x_new, method = "linear")$y
  }
  
  return(interpolated)
}

# Function to prepare data for WICMAD
prepare_wicmad_data <- function(series_list, labels) {
  cat("Preparing data for WICMAD...\n")
  
  n_series <- length(series_list)
  
  # Find the highest power of 2 supported by the data
  max_dim <- max(sapply(series_list, nrow))
  p_dim <- 2^floor(log2(max_dim))
  cat("Using p parameter (highest power of 2):", p_dim, "\n")
  
  # Interpolate each series to p dimensions
  interpolated_series <- list()
  for (i in seq_len(n_series)) {
    interpolated_series[[i]] <- interpolate_series(series_list[[i]], p_dim)
  }
  
  # Create time coordinates
  t <- seq(0, 1, length.out = p_dim)
  
  # Convert to list format for WICMAD (each element is a PxM matrix)
  # P=p_dim (time points), M=6 (sensors)
  Y <- interpolated_series
  
  return(list(Y = Y, t = t, labels = labels))
}

# Function to calculate clustering metrics
calculate_clustering_metrics <- function(true_labels, pred_labels) {
  true <- factor(true_labels, levels = sort(unique(c(true_labels, pred_labels))))
  pred <- factor(pred_labels, levels = sort(unique(c(true_labels, pred_labels))))
  
  cm <- table(true, pred)
  cat("Confusion Matrix:\n")
  print(cm)
  
  # Calculate metrics
  accuracy <- sum(diag(cm)) / sum(cm)
  precision <- diag(cm) / rowSums(cm)
  recall <- diag(cm) / colSums(cm)
  f1 <- 2 * (precision * recall) / (precision + recall)
  
  # Adjusted Rand Index
  ari <- adjustedRandIndex(true_labels, pred_labels)
  
  cat("Accuracy:", round(accuracy, 3), "\n")
  cat("Precision:", round(precision, 3), "\n")
  cat("Recall:", round(recall, 3), "\n")
  cat("F1-Score:", round(f1, 3), "\n")
  cat("Adjusted Rand Index:", round(ari, 3), "\n")
  
  return(list(
    Confusion_Matrix = cm,
    Accuracy = accuracy,
    Precision = precision,
    Recall = recall,
    F1_Score = f1,
    ARI = ari
  ))
}

# Function to plot overlapped data
plot_overlapped_data <- function(series_list, labels, title = "UCI HAR - Before Clustering", max_series = NULL) {
  cat("Creating overlapped data plot...\n")
  
  # Use all series if max_series is NULL
  if (is.null(max_series)) {
    series_indices <- seq_len(length(series_list))
  } else {
    n_series <- min(max_series, length(series_list))
    series_indices <- sample(seq_len(length(series_list)), n_series)
  }
  
  # Create plotting data
  plot_data <- data.frame()
  
  for (i in series_indices) {
    series_data <- series_list[[i]]
    time_points <- seq(0, 1, length.out = nrow(series_data))
    
    for (j in seq_len(ncol(series_data))) {
      series_df <- data.frame(
        time = time_points,
        value = series_data[, j],
        dimension = paste("Sensor", j),
        dimension_order = j,  # Add numeric ordering
        series_id = i,
        label = labels[i],
        label_name = ifelse(labels[i] == 0, "Normal", "Anomaly"),
        row.names = NULL
      )
      
      plot_data <- rbind(plot_data, series_df)
    }
  }
  
  # Create ordered factor for proper sensor ordering
  plot_data$dimension <- factor(plot_data$dimension, 
                                 levels = paste("Sensor", sort(unique(plot_data$dimension_order))))
  
  p <- ggplot(plot_data, aes(x = .data$time, y = .data$value, color = .data$label_name, group = .data$series_id)) +
    geom_line(alpha = 0.6, linewidth = 0.5) +
    facet_wrap(~ .data$dimension, scales = "free_y", ncol = 3) +
    labs(title = title,
         x = "Time", y = "Sensor Value", color = "Class") +
    scale_color_manual(values = c("Normal" = "blue", "Anomaly" = "red")) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(p)
}

# Function to plot clustering results
plot_clustering_results <- function(series_list, true_labels, cluster_assignments, title = "UCI HAR - After Clustering", max_series = NULL) {
  cat("Creating clustering results plot...\n")
  
  # Use all series if max_series is NULL
  if (is.null(max_series)) {
    series_indices <- seq_len(length(series_list))
  } else {
    n_series <- min(max_series, length(series_list))
    series_indices <- sample(seq_len(length(series_list)), n_series)
  }
  
  # Create plotting data
  plot_data <- data.frame()
  
  for (i in series_indices) {
    series_data <- series_list[[i]]
    time_points <- seq(0, 1, length.out = nrow(series_data))
    
    for (j in seq_len(ncol(series_data))) {
      series_df <- data.frame(
        time = time_points,
        value = series_data[, j],
        dimension = paste("Sensor", j),
        dimension_order = j,  # Add numeric ordering
        series_id = i,
        true_label = true_labels[i],
        cluster = cluster_assignments[i],
        cluster_name = ifelse(cluster_assignments[i] == 0, "Normal", "Anomaly"),
        row.names = NULL
      )
      
      plot_data <- rbind(plot_data, series_df)
    }
  }
  
  # Create ordered factor for proper sensor ordering
  plot_data$dimension <- factor(plot_data$dimension, 
                                 levels = paste("Sensor", sort(unique(plot_data$dimension_order))))
  
  p <- ggplot(plot_data, aes(x = .data$time, y = .data$value, color = .data$cluster_name, group = .data$series_id)) +
    geom_line(alpha = 0.6, linewidth = 0.5) +
    facet_wrap(~ .data$dimension, scales = "free_y", ncol = 3) +
    labs(title = title,
         x = "Time", y = "Sensor Value", color = "Class") +
    scale_color_manual(values = c("Normal" = "blue", "Anomaly" = "red")) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(p)
}

# Function to run WICMAD analysis
run_wicmad_analysis <- function(prepare_func, analysis_name, series_data, labels) {
  cat(paste("\n=== Running WICMAD analysis on", analysis_name, "===\n"))
  
  # Prepare data
  wicmad_data <- prepare_func(series_data, labels)
  
  # Run WICMAD with test parameters
  cat("Running WICMAD clustering...\n")
  cat("Parameters: n_iter=10, burn=5\n")
  
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
  
  # Map WICMAD clusters to binary classification (normal vs anomaly)
  # Use pinned cluster (0) as normal, all others as anomaly
  binary_cluster_assignments <- ifelse(cluster_assignments == 0, 0, 1)
  
  cat("WICMAD found", length(unique(cluster_assignments)), "clusters\n")
  cat("Pinned cluster (0) -> normal, all others -> anomaly\n")
  
  # Calculate metrics
  metrics <- calculate_clustering_metrics(labels, binary_cluster_assignments)
  
  return(list(
    wicmad_result = wicmad_result,
    cluster_assignments = binary_cluster_assignments,
    metrics = metrics
  ))
}

# Main function
main <- function() {
  cat("=== UCI HAR Dataset Analysis with WICMAD ===\n")
  
  # Create output directory
  output_dir <- "../../plots/ucihar"
  dir.create(output_dir, recursive = TRUE)
  
  # Load data
  cat("\n1. Loading UCI HAR dataset with 10% anomaly class...\n")
  data <- load_ucihar_data(data_dir = "../../data")
  
  # Create original data visualization
  cat("\n2. Creating original data visualization...\n")
  
  # Plot 500 normal curves and all anomaly curves for faster visualization
  normal_indices <- which(data$train_labels == 0)
  anomaly_indices <- which(data$train_labels == 1)
  
  # Sample 500 normal curves (or all if less than 500)
  n_normal_to_plot <- min(500, length(normal_indices))
  selected_normal <- sample(normal_indices, n_normal_to_plot)
  
  # Use all anomaly curves
  plot_indices <- c(selected_normal, anomaly_indices)
  plot_series <- data$train_series[plot_indices]
  plot_labels <- data$train_labels[plot_indices]
  
  cat("Plotting", length(selected_normal), "normal curves and", length(anomaly_indices), "anomaly curves (", length(plot_series), "total)\n")
  
  original_plot <- plot_overlapped_data(plot_series, plot_labels, 
                                       "UCI HAR - Before Clustering")
  
  # Save original data plot
  pdf("../../plots/ucihar/ucihar_original_data.pdf", width = 12, height = 8)
  print(original_plot)
  dev.off()
  cat("Original data plot saved to ../../plots/ucihar/ucihar_original_data.pdf\n")
  
  # Run WICMAD analysis on raw data
  cat("\n3. Analyzing raw sensor data...\n")
  raw_results <- run_wicmad_analysis(
    function(series, labels) prepare_wicmad_data(series, labels), 
    "Raw Sensor Data", 
    data$train_series, 
    data$train_labels
  )
  
  # Plot clustering results (using same subset as original plot: 500 normal + all anomalies)
  plot_cluster_series <- data$train_series[plot_indices]
  plot_cluster_labels <- data$train_labels[plot_indices]
  plot_cluster_assignments <- raw_results$cluster_assignments[plot_indices]
  
  clustering_plot <- plot_clustering_results(
    plot_cluster_series, 
    plot_cluster_labels, 
    plot_cluster_assignments,
    "UCI HAR - After Clustering"
  )
  
  # Save clustering results plot
  pdf("../../plots/ucihar/ucihar_clustering_results.pdf", width = 12, height = 8)
  print(clustering_plot)
  dev.off()
  cat("Clustering results plot saved to ../../plots/ucihar/ucihar_clustering_results.pdf\n")
  
  # Print final results
  cat("\n=== Final Results ===\n")
  cat("Raw Sensor Data Analysis:\n")
  cat("Accuracy:", round(raw_results$metrics$Accuracy, 3), "\n")
  cat("Adjusted Rand Index:", round(raw_results$metrics$ARI, 3), "\n")
  
  cat("\nAnalysis complete! Plots saved to ../../plots/ucihar/\n")
}

# Run the analysis
main()