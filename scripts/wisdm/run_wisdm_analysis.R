#!/usr/bin/env Rscript

# WISDM Dataset Analysis with WICMAD
# This script loads the WISDM dataset, runs WICMAD clustering, and visualizes results

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(cluster)
library(mclust)  # For adjustedRandIndex

# Load WICMAD package
library(devtools)
devtools::load_all()

# Function to load WISDM dataset and create imbalanced data
load_wisdm_data <- function(data_dir) {
  cat("Loading WISDM dataset...\n")
  
  # Load activity key
  activity_key <- read.table(file.path(data_dir, "wisdm-dataset", "wisdm-dataset", "activity_key.txt"), 
                             sep = "=", col.names = c("activity", "code"), strip.white = TRUE)
  
  # Load raw accelerometer data from phone
  cat("Loading accelerometer data...\n")
  accel_files <- list.files(file.path(data_dir, "wisdm-dataset", "wisdm-dataset", "raw", "phone", "accel"), 
                            pattern = "*.txt", full.names = TRUE)
  
  # Load raw gyroscope data from phone
  cat("Loading gyroscope data...\n")
  gyro_files <- list.files(file.path(data_dir, "wisdm-dataset", "wisdm-dataset", "raw", "phone", "gyro"), 
                          pattern = "*.txt", full.names = TRUE)
  
  # Function to read sensor data file
  read_sensor_file <- function(file_path) {
    # Read the file line by line and clean the semicolon
    lines <- readLines(file_path)
    # Remove semicolon from the end of each line
    lines <- gsub(";$", "", lines)
    # Convert to data frame
    data <- read.table(text = lines, sep = ",", header = FALSE, 
                      col.names = c("user", "activity", "timestamp", "x", "y", "z"))
    return(data)
  }
  
  # Load all accelerometer data
  cat("Processing accelerometer data...\n")
  accel_data <- list()
  for (file in accel_files[1:min(10, length(accel_files))]) {  # Limit to first 10 files for speed
    accel_data[[length(accel_data) + 1]] <- read_sensor_file(file)
  }
  accel_combined <- do.call(rbind, accel_data)
  
  # Load all gyroscope data
  cat("Processing gyroscope data...\n")
  gyro_data <- list()
  for (file in gyro_files[1:min(10, length(gyro_files))]) {  # Limit to first 10 files for speed
    gyro_data[[length(gyro_data) + 1]] <- read_sensor_file(file)
  }
  gyro_combined <- do.call(rbind, gyro_data)
  
  # Create time series segments
  cat("Creating time series segments...\n")
  segment_length <- 60  # 60 time points per segment
  time_series_data <- list()
  labels <- c()
  
  # Get unique activities
  unique_activities <- unique(accel_combined$activity)
  cat("Found activities:", paste(unique_activities, collapse = ", "), "\n")
  
  # Create segments for each activity
  for (activity in unique_activities) {
    # Get data for this activity
    accel_activity <- accel_combined[accel_combined$activity == activity, ]
    gyro_activity <- gyro_combined[gyro_combined$activity == activity, ]
    
    # Create segments
    n_segments <- min(50, floor(nrow(accel_activity) / segment_length))  # Limit segments per activity
    
    for (i in seq_len(n_segments)) {
      start_idx <- (i - 1) * segment_length + 1
      end_idx <- min(start_idx + segment_length - 1, nrow(accel_activity))
      
      if (end_idx - start_idx + 1 >= segment_length) {
        # Create multivariate time series (6 sensors: 3 accel + 3 gyro)
        series_matrix <- matrix(0, nrow = segment_length, ncol = 6)
        
        # Fill accelerometer data
        series_matrix[, 1] <- accel_activity$x[start_idx:end_idx]
        series_matrix[, 2] <- accel_activity$y[start_idx:end_idx]
        series_matrix[, 3] <- accel_activity$z[start_idx:end_idx]
        
        # Fill gyroscope data (if available)
        if (nrow(gyro_activity) >= end_idx) {
          series_matrix[, 4] <- gyro_activity$x[start_idx:end_idx]
          series_matrix[, 5] <- gyro_activity$y[start_idx:end_idx]
          series_matrix[, 6] <- gyro_activity$z[start_idx:end_idx]
        }
        
        time_series_data[[length(time_series_data) + 1]] <- series_matrix
        labels <- c(labels, activity)
      }
    }
  }
  
  cat("Dataset dimensions:\n")
  cat("Number of samples:", length(time_series_data), "\n")
  cat("Time points per sample:", nrow(time_series_data[[1]]), "\n")
  cat("Number of sensors:", ncol(time_series_data[[1]]), "\n")
  
  # Find the majority class to use as normal
  class_counts <- table(labels)
  majority_class <- names(class_counts)[which.max(class_counts)]
  minority_classes <- names(class_counts)[class_counts != max(class_counts)]
  
  cat("Class distribution:\n")
  print(class_counts)
  cat("Majority class (normal):", majority_class, "\n")
  cat("Minority classes (anomalies):", paste(minority_classes, collapse = ", "), "\n")
  
  # Convert to anomaly detection format (majority class = normal, others = anomaly)
  anomaly_labels <- ifelse(labels == majority_class, 0, 1)
  
  # Create imbalanced dataset with 5% anomalies using all available data
  normal_indices <- which(anomaly_labels == 0)
  anomaly_indices <- which(anomaly_labels == 1)
  
  # Calculate how many anomalies we need for 5% of total
  total_samples <- length(anomaly_labels)
  n_anomalies_needed <- max(1, round(total_samples * 0.05))
  n_normal_needed <- total_samples - n_anomalies_needed
  
  cat("Creating imbalanced dataset with 5% anomalies:\n")
  cat("Total samples available:", total_samples, "\n")
  cat("Normal samples needed:", n_normal_needed, "\n")
  cat("Anomaly samples needed:", n_anomalies_needed, "\n")
  
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
  imbalanced_series <- time_series_data[selected_indices]
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
plot_overlapped_data <- function(series_list, labels, title = "WISDM - Before Clustering", max_series = NULL) {
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
plot_clustering_results <- function(series_list, true_labels, cluster_assignments, title = "WISDM - After Clustering", max_series = NULL) {
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
  cat("=== WISDM Dataset Analysis with WICMAD ===\n")
  
  # Create output directory
  output_dir <- "../plots/wisdm"
  dir.create(output_dir, recursive = TRUE)
  
  # Load data
  cat("\n1. Loading WISDM dataset with 5% anomaly class...\n")
  data <- load_wisdm_data(data_dir = "../data")
  
  # Create original data visualization
  cat("\n2. Creating original data visualization...\n")
  
  # Use all samples for plotting
  plot_indices <- seq_len(length(data$train_series))
  plot_series <- data$train_series[plot_indices]
  plot_labels <- data$train_labels[plot_indices]
  
  cat("Plotting all", length(plot_series), "curves\n")
  
  original_plot <- plot_overlapped_data(plot_series, plot_labels, 
                                       "WISDM - Before Clustering")
  
  # Save original data plot
  pdf("../plots/wisdm/wisdm_original_data.pdf", width = 12, height = 8)
  print(original_plot)
  dev.off()
  cat("Original data plot saved to ../plots/wisdm/wisdm_original_data.pdf\n")
  
  # Run WICMAD analysis on raw data
  cat("\n3. Analyzing raw sensor data...\n")
  raw_results <- run_wicmad_analysis(
    function(series, labels) prepare_wicmad_data(series, labels), 
    "Raw Sensor Data", 
    data$train_series, 
    data$train_labels
  )
  
  # Plot clustering results
  plot_cluster_series <- data$train_series[plot_indices]
  plot_cluster_labels <- data$train_labels[plot_indices]
  plot_cluster_assignments <- raw_results$cluster_assignments[plot_indices]
  
  clustering_plot <- plot_clustering_results(
    plot_cluster_series, 
    plot_cluster_labels, 
    plot_cluster_assignments,
    "WISDM - After Clustering"
  )
  
  # Save clustering results plot
  pdf("../plots/wisdm/wisdm_clustering_results.pdf", width = 12, height = 8)
  print(clustering_plot)
  dev.off()
  cat("Clustering results plot saved to ../plots/wisdm/wisdm_clustering_results.pdf\n")
  
  # Print final results
  cat("\n=== Final Results ===\n")
  cat("Raw Sensor Data Analysis:\n")
  cat("Accuracy:", round(raw_results$metrics$Accuracy, 3), "\n")
  cat("Adjusted Rand Index:", round(raw_results$metrics$ARI, 3), "\n")
  
  cat("\nAnalysis complete! Plots saved to ../plots/wisdm/\n")
}

# Run the analysis
main()