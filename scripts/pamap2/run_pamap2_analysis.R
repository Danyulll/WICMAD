#!/usr/bin/env Rscript

# PAMAP2 Dataset Analysis with WICMAD
# This script loads the PAMAP2 dataset, runs WICMAD clustering, and visualizes results

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(cluster)
library(mclust)  # For adjustedRandIndex

# Load WICMAD package
library(devtools)
devtools::load_all()

# Function to load PAMAP2 dataset and create imbalanced data
load_pamap2_data <- function(data_dir) {
  cat("Loading PAMAP2 dataset...\n")
  
  # PAMAP2 has 18 activities, we'll use the most common as normal
  # Activities: 1=Lying, 2=Sitting, 3=Standing, 4=Walking, 5=Running, 6=Cycling, 
  # 7=Nordic walking, 8=Ascending stairs, 9=Descending stairs, 10=Vacuum cleaning,
  # 11=Ironing, 12=Rope jumping, 13=Watching TV, 14=Computer work, 15=Car driving,
  # 16=Standing up from lying, 17=Lying down from standing, 18=Transition from sitting to standing
  
  # Load real PAMAP2 data from Protocol directory
  protocol_dir <- file.path(data_dir, "PAMAP2_Dataset", "PAMAP2_Dataset", "Protocol")
  
  if (!dir.exists(protocol_dir)) {
    stop("PAMAP2 Protocol directory not found: ", protocol_dir)
  }
  
  # Get all subject files
  subject_files <- list.files(protocol_dir, pattern = "subject.*\\.dat", full.names = TRUE)
  cat("Found", length(subject_files), "subject files\n")
  
  # Function to read a subject file
  read_subject_file <- function(file_path) {
    cat("Reading", basename(file_path), "...\n")
    
    # Read the data file
    data <- read.table(file_path, header = FALSE, na.strings = "NaN")
    
    # PAMAP2 data format: timestamp, activity_id, heart_rate, temperature, 
    # then 3 IMU sensors (hand, chest, ankle) with 3 acc + 3 gyro + 3 mag each
    # Total: 1 + 1 + 1 + 1 + 3*(3+3+3) = 1 + 1 + 1 + 1 + 27 = 31 columns
    
    # Extract relevant columns (skip timestamp, use activity_id, heart_rate, and IMU data)
    # Keep: activity_id (col 2), heart_rate (col 3), and IMU data (cols 5-31)
    activity_col <- 2
    heart_rate_col <- 3
    imu_start_col <- 5
    imu_end_col <- 31
    
    # Extract data
    activities <- data[, activity_col]
    heart_rates <- data[, heart_rate_col]
    imu_data <- data[, imu_start_col:imu_end_col]
    
    # Remove rows with NaN values
    valid_rows <- !is.na(activities) & !is.na(heart_rates) & 
                  rowSums(is.na(imu_data)) == 0
    
    if (sum(valid_rows) == 0) {
      cat("Warning: No valid data in", basename(file_path), "\n")
      return(NULL)
    }
    
    activities <- activities[valid_rows]
    heart_rates <- heart_rates[valid_rows]
    imu_data <- imu_data[valid_rows, ]
    
    return(list(
      activities = activities,
      heart_rates = heart_rates,
      imu_data = imu_data
    ))
  }
  
  # Load all subject data
  all_data <- list()
  for (file in subject_files) {
    subject_data <- read_subject_file(file)
    if (!is.null(subject_data)) {
      all_data[[length(all_data) + 1]] <- subject_data
    }
  }
  
  if (length(all_data) == 0) {
    stop("No valid PAMAP2 data found")
  }
  
  # Combine all subject data
  all_activities <- unlist(lapply(all_data, function(x) x$activities))
  all_heart_rates <- unlist(lapply(all_data, function(x) x$heart_rates))
  all_imu_data <- do.call(rbind, lapply(all_data, function(x) x$imu_data))
  
  cat("Combined dataset:\n")
  cat("Total samples:", length(all_activities), "\n")
  cat("IMU data dimensions:", nrow(all_imu_data), "x", ncol(all_imu_data), "\n")
  
  # Create time series segments from the data
  # Use sliding window approach to create segments
  segment_length <- 100  # 100 time points per segment
  time_series_data <- list()
  segment_labels <- c()
  
  # Create segments for each activity
  for (activity in unique(all_activities)) {
    activity_indices <- which(all_activities == activity)
    
    if (length(activity_indices) < segment_length) {
      next  # Skip activities with insufficient data
    }
    
    # Create overlapping segments
    n_segments <- min(50, floor(length(activity_indices) / segment_length))  # Limit segments per activity
    
    for (i in seq_len(n_segments)) {
      start_idx <- (i - 1) * (segment_length / 2) + 1  # 50% overlap
      end_idx <- min(start_idx + segment_length - 1, length(activity_indices))
      
      if (end_idx - start_idx + 1 >= segment_length) {
        # Extract segment data
        segment_indices <- activity_indices[start_idx:end_idx]
        
        # Create multivariate time series (heart_rate + IMU data)
        # Combine heart rate and IMU data
        segment_data <- cbind(
          all_heart_rates[segment_indices],
          all_imu_data[segment_indices, ]
        )
        
        time_series_data[[length(time_series_data) + 1]] <- segment_data
        segment_labels <- c(segment_labels, activity)
      }
    }
  }
  
  cat("Created", length(time_series_data), "time series segments\n")
  
  # Find the majority class to use as normal
  class_counts <- table(segment_labels)
  majority_class <- as.numeric(names(class_counts)[which.max(class_counts)])
  minority_classes <- as.numeric(names(class_counts)[class_counts != max(class_counts)])
  
  cat("Class distribution:\n")
  print(class_counts)
  cat("Majority class (normal):", majority_class, "\n")
  cat("Minority classes (anomalies):", paste(minority_classes, collapse = ", "), "\n")
  
  # Convert segment labels to anomaly format
  segment_anomaly_labels <- ifelse(segment_labels == majority_class, 0, 1)
  
  # Create imbalanced dataset with approximately 5% anomalies
  normal_indices <- which(segment_anomaly_labels == 0)
  anomaly_indices <- which(segment_anomaly_labels == 1)
  
  # Use all available normal samples
  selected_normal <- normal_indices
  n_normal_used <- length(selected_normal)
  
  # Calculate how many anomalies we need for 5% of the total dataset
  # Total dataset will be: n_normal_used + n_anomalies_needed
  # We want: n_anomalies_needed / (n_normal_used + n_anomalies_needed) = 0.05
  # Solving: n_anomalies_needed = 0.05 * (n_normal_used + n_anomalies_needed)
  # n_anomalies_needed = 0.05 * n_normal_used / (1 - 0.05) = 0.05 * n_normal_used / 0.95
  n_anomalies_needed <- max(1, round(0.05 * n_normal_used / 0.95))
  
  # Ensure we don't exceed available anomaly samples
  n_anomalies_needed <- min(n_anomalies_needed, length(anomaly_indices))
  
  # Sample the required number of anomaly samples
  selected_anomaly <- sample(anomaly_indices, n_anomalies_needed)
  
  # Combine selected indices
  selected_indices <- c(selected_normal, selected_anomaly)
  
  # Create imbalanced dataset
  imbalanced_series <- time_series_data[selected_indices]
  imbalanced_labels <- segment_anomaly_labels[selected_indices]
  
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
  # P=p_dim (time points), M=9 (sensors)
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
plot_overlapped_data <- function(series_list, labels, title = "PAMAP2 - Before Clustering", max_series = NULL) {
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
plot_clustering_results <- function(series_list, true_labels, cluster_assignments, title = "PAMAP2 - After Clustering", max_series = NULL) {
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
  cat("Parameters: n_iter=8000, burn=3000, thin=5\n")
  
  # Run WICMAD
  wicmad_result <- wicmad(
    Y = wicmad_data$Y,
    t = wicmad_data$t,
    n_iter = 8000,
    burn = 3000,
    thin = 5
  )
  
  # Extract cluster assignments
  cluster_assignments <- wicmad_result$Z[nrow(wicmad_result$Z), ]
  
  # Map WICMAD clusters to binary classification (normal vs anomaly)
  # Use the most frequent cluster as normal, all others as anomaly
  cluster_counts <- table(cluster_assignments)
  normal_cluster <- as.numeric(names(cluster_counts)[which.max(cluster_counts)])
  binary_cluster_assignments <- ifelse(cluster_assignments == normal_cluster, 0, 1)
  
  cat("WICMAD found", length(unique(cluster_assignments)), "clusters\n")
  cat("Most frequent cluster (", normal_cluster, ") -> normal, all others -> anomaly\n")
  
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
  cat("=== PAMAP2 Dataset Analysis with WICMAD ===\n")
  
  # Create output directory
  output_dir <- "../../plots/pamap2"
  dir.create(output_dir, recursive = TRUE)
  
  # Load data
  cat("\n1. Loading PAMAP2 dataset with 5% anomaly class...\n")
  data <- load_pamap2_data(data_dir = "../../data")
  
  # Create original data visualization
  cat("\n2. Creating original data visualization...\n")
  
  # Use a subset for plotting (first 50 series)
  plot_indices <- seq_len(min(50, length(data$train_series)))
  plot_series <- data$train_series[plot_indices]
  plot_labels <- data$train_labels[plot_indices]
  
  cat("Plotting subsample: 50 curves\n")
  
  original_plot <- plot_overlapped_data(plot_series, plot_labels, 
                                       "PAMAP2 - Before Clustering")
  
  # Save original data plot
  pdf("../../plots/pamap2/pamap2_original_data.pdf", width = 12, height = 8)
  print(original_plot)
  dev.off()
  cat("Original data plot saved to ../../plots/pamap2/pamap2_original_data.pdf\n")
  
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
    "PAMAP2 - After Clustering"
  )
  
  # Save clustering results plot
  pdf("../../plots/pamap2/pamap2_clustering_results.pdf", width = 12, height = 8)
  print(clustering_plot)
  dev.off()
  cat("Clustering results plot saved to ../../plots/pamap2/pamap2_clustering_results.pdf\n")
  
  # Print final results
  cat("\n=== Final Results ===\n")
  cat("Raw Sensor Data Analysis:\n")
  cat("Accuracy:", round(raw_results$metrics$Accuracy, 3), "\n")
  cat("Adjusted Rand Index:", round(raw_results$metrics$ARI, 3), "\n")
  
  cat("\nAnalysis complete! Plots saved to ../../plots/pamap2/\n")
}

# Run the analysis
main()