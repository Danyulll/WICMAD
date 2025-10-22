#!/usr/bin/env Rscript

# BasicMotions Dataset Analysis with WICMAD
# This script loads the BasicMotions dataset, runs WICMAD clustering, and visualizes results

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(cluster)
# library(caret)
library(mclust)  # For adjustedRandIndex

# Load WICMAD package
library(devtools)
devtools::install()

# Function to load BasicMotions dataset with imbalanced sampling
load_basicmotions_data <- function(data_dir, reveal_ratio = 0.15) {
  cat("Loading BasicMotions dataset...\n")
  
  # Load training data
  train_file <- file.path(data_dir, "BasicMotions_TRAIN.ts")
  test_file <- file.path(data_dir, "BasicMotions_TEST.ts")
  
  if (!file.exists(train_file)) {
    stop("Training file not found: ", train_file)
  }
  
  # Read the .ts file
  train_data <- readLines(train_file)
  
  # Find the @data line
  data_start <- which(grepl("^@data", train_data))
  if (length(data_start) == 0) {
    stop("No @data section found in file")
  }
  
  # Extract data lines
  data_lines <- train_data[(data_start + 1):length(train_data)]
  data_lines <- data_lines[data_lines != ""]  # Remove empty lines
  
  # Parse the data
  n_series <- length(data_lines)
  cat("Number of series:", n_series, "\n")
  
  # Parse each series
  series_list <- list()
  labels <- character(n_series)
  
  for (i in seq_len(n_series)) {
    line <- data_lines[i]
    
    # Split by colon to separate data from label
    parts <- strsplit(line, ":")[[1]]
    if (length(parts) != 7) {
      stop("Invalid data format in line ", i, " - expected 7 parts, got ", length(parts))
    }
    
    # Parse the data parts (6 dimensions)
    data_values <- numeric()
    for (j in 1:6) {
      dim_values <- as.numeric(strsplit(parts[j], ",")[[1]])
      data_values <- c(data_values, dim_values)
    }
    
    # Parse the label part (string label)
    labels[i] <- parts[7]
    
    # Reshape data into 6 columns (dimensions)
    n_timepoints <- length(data_values) / 6
    series_matrix <- matrix(data_values, nrow = n_timepoints, ncol = 6, byrow = FALSE)
    
    series_list[[i]] <- series_matrix
  }
  
  # Convert string labels to numeric
  unique_labels <- unique(labels)
  label_mapping <- setNames(seq_len(length(unique_labels)) - 1, unique_labels)
  labels <- label_mapping[labels]
  
  # Load test data if available
  test_series_list <- NULL
  test_labels <- NULL
  if (file.exists(test_file)) {
    cat("Loading test data...\n")
    test_data <- readLines(test_file)
    test_data_start <- which(grepl("^@data", test_data))
    test_data_lines <- test_data[(test_data_start + 1):length(test_data)]
    test_data_lines <- test_data_lines[test_data_lines != ""]
    
    n_test_series <- length(test_data_lines)
    test_series_list <- list()
    test_labels <- integer(n_test_series)
    
    for (i in seq_len(n_test_series)) {
      line <- test_data_lines[i]
      parts <- strsplit(line, ":")[[1]]
      
      # Parse the data parts (6 dimensions)
      data_values <- numeric()
      for (j in 1:6) {
        dim_values <- as.numeric(strsplit(parts[j], ",")[[1]])
        data_values <- c(data_values, dim_values)
      }
      
      # Parse the label part (string label)
      test_labels[i] <- parts[7]
      
      n_timepoints <- length(data_values) / 6
      series_matrix <- matrix(data_values, nrow = n_timepoints, ncol = 6, byrow = FALSE)
      
      test_series_list[[i]] <- series_matrix
    }
    
    # Convert test labels to numeric using same mapping
    test_labels <- label_mapping[test_labels]
  }
  
  # Combine train and test data
  all_series <- c(series_list, test_series_list)
  all_labels <- c(labels, test_labels)
  
  cat("Combined dataset:\n")
  cat("Train series:", length(series_list), "\n")
  cat("Test series:", length(test_series_list), "\n")
  cat("Total series:", length(all_series), "\n")
  
  # Create imbalanced dataset: Walking (normal) + 1 from each other class
  cat("Creating imbalanced dataset...\n")
  cat("Walking as normal class, 1 observation from each other class\n")
  
  # Find Walking class (assuming it's one of the labels)
  walking_indices <- which(all_labels == which(unique_labels == "Walking") - 1)
  other_indices <- which(all_labels != which(unique_labels == "Walking") - 1)
  
  cat("Walking indices found:", length(walking_indices), "\n")
  cat("Other indices found:", length(other_indices), "\n")
  
  # Sample 1 from each other class
  other_classes <- unique(all_labels[other_indices])
  selected_other <- integer()
  for (class in other_classes) {
    class_indices <- which(all_labels == class)
    selected_other <- c(selected_other, sample(class_indices, 1))
  }
  
  # Use all Walking samples as normal class (no sampling)
  revealed_walking <- walking_indices
  
  # Combine selected samples
  selected_indices <- c(revealed_walking, selected_other)
  
  imbalanced_series <- all_series[selected_indices]
  imbalanced_labels <- all_labels[selected_indices]
  
  # Convert to anomaly detection format (0 = normal/Walking, 1 = anomaly/other)
  anomaly_labels <- ifelse(imbalanced_labels == which(unique_labels == "Walking") - 1, 0, 1)
  
  cat("Final imbalanced dataset:\n")
  cat("Total samples:", length(imbalanced_series), "\n")
  cat("Normal (Walking):", sum(anomaly_labels == 0), "\n")
  cat("Anomaly (Other):", sum(anomaly_labels == 1), "\n")
  cat("Anomaly percentage:", round(mean(anomaly_labels == 1) * 100, 1), "%\n")
  
  return(list(
    train_series = imbalanced_series,
    train_labels = anomaly_labels,
    test_series = NULL,
    test_labels = NULL
  ))
}

# Function to interpolate time series to target dimension
interpolate_series <- function(series, target_dim = 64) {
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
  
  # Interpolate each series to 64 dimensions (nearest power of 2 to 100)
  interpolated_series <- list()
  for (i in seq_len(n_series)) {
    interpolated_series[[i]] <- interpolate_series(series_list[[i]], 64)
  }
  
  # Create time coordinates
  t <- seq(0, 1, length.out = 64)
  
  # Convert to list format for WICMAD (each element is a PxM matrix)
  # P=64 (time points), M=6 (dimensions)
  Y <- interpolated_series
  
  return(list(Y = Y, t = t, labels = labels))
}

# Function to calculate clustering metrics
calculate_clustering_metrics <- function(true_labels, pred_labels) {
  # Convert to factors with same levels
  true <- factor(true_labels, levels = sort(unique(c(true_labels, pred_labels))))
  pred <- factor(pred_labels, levels = sort(unique(c(true_labels, pred_labels))))
  
  # Create confusion matrix
  cm <- table(true, pred)
  cat("Confusion Matrix:\n")
  print(cm)
  
  # Calculate metrics for each class
  classes <- levels(true)
  metrics <- list()
  
  for (class in classes) {
    tp <- sum(true == class & pred == class)
    fp <- sum(true != class & pred == class)
    fn <- sum(true == class & pred != class)
    tn <- sum(true != class & pred != class)
    
    # Calculate metrics
    precision <- if (tp + fp == 0) 0 else tp / (tp + fp)
    recall <- if (tp + fn == 0) 0 else tp / (tp + fn)
    specificity <- if (tn + fp == 0) 0 else tn / (tn + fp)
    f1 <- if (precision + recall == 0) 0 else 2 * precision * recall / (precision + recall)
    accuracy <- (tp + tn) / (tp + tn + fp + fn)
    
    metrics[[paste0("Class_", class)]] <- list(
      Precision = precision,
      Recall = recall,
      Specificity = specificity,
      F1 = f1,
      Accuracy = accuracy
    )
  }
  
  # Calculate macro-averaged metrics
  macro_precision <- mean(sapply(metrics, function(x) x$Precision), na.rm = TRUE)
  macro_recall <- mean(sapply(metrics, function(x) x$Recall), na.rm = TRUE)
  macro_specificity <- mean(sapply(metrics, function(x) x$Specificity), na.rm = TRUE)
  macro_f1 <- mean(sapply(metrics, function(x) x$F1), na.rm = TRUE)
  overall_accuracy <- mean(sapply(metrics, function(x) x$Accuracy), na.rm = TRUE)
  
  # Adjusted Rand Index
  ari <- adjustedRandIndex(true_labels, pred_labels)
  
  return(list(
    Confusion_Matrix = cm,
    Macro_Precision = macro_precision,
    Macro_Recall = macro_recall,
    Macro_Specificity = macro_specificity,
    Macro_F1 = macro_f1,
    Overall_Accuracy = overall_accuracy,
    ARI = ari,
    Class_Metrics = metrics
  ))
}

# Function to plot overlapped data
plot_overlapped_data <- function(series_list, labels, title = "BasicMotions Dataset", max_series = NULL) {
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
        dimension = paste("Dim", j),
        series_id = i,
        label = labels[i],
        label_name = paste("Class", labels[i]),
        row.names = NULL
      )
      
      plot_data <- rbind(plot_data, series_df)
    }
  }
  
  # Create overlapped plot
  p <- ggplot(plot_data, aes(x = .data$time, y = .data$value, color = .data$label_name, group = .data$series_id)) +
    geom_line(alpha = 0.6, linewidth = 0.5) +
    facet_wrap(~ .data$dimension, scales = "free_y", ncol = 2) +
    labs(title = title,
         x = "Time", y = "Value", color = "Class") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(p)
}

# Function to plot clustering results
plot_clustering_results <- function(series_list, true_labels, cluster_assignments, title = "BasicMotions - Clustering Results", max_series = NULL) {
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
        dimension = paste("Dim", j),
        series_id = i,
        true_label = true_labels[i],
        cluster = cluster_assignments[i],
        true_label_name = paste("Class", true_labels[i]),
        row.names = NULL
      )
      
      plot_data <- rbind(plot_data, series_df)
    }
  }
  
  # Create overlapped plot
  p <- ggplot(plot_data, aes(x = .data$time, y = .data$value, color = .data$true_label_name, group = .data$series_id)) +
    geom_line(alpha = 0.6, linewidth = 0.5) +
    facet_wrap(~ .data$dimension, scales = "free_y", ncol = 2) +
    labs(title = title,
         x = "Time", y = "Value", color = "True Class") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(p)
}

# Main analysis
main <- function() {
  cat("=== BasicMotions Dataset Analysis with WICMAD ===\n\n")
  
  # Set data directory
  data_dir <- "data/BasicMotions"
  
  # Load data
  cat("1. Loading BasicMotions dataset...\n")
  data <- load_basicmotions_data(data_dir, reveal_ratio = 0.15)
  
  # Print data summary
  cat("Training data: ", length(data$train_series), " series\n")
  cat("Label distribution:\n")
  print(table(data$train_labels))
  cat("First few labels:", head(data$train_labels, 10), "\n")
  
  # Plot original data (overlapped)
  cat("\n2. Creating original data visualization...\n")
  original_plot <- plot_overlapped_data(data$train_series, data$train_labels, 
                                       "BasicMotions Dataset - Original Multivariate Time Series")
  
  # Save original data plot
  pdf("basicmotions_original_data.pdf", width = 12, height = 10)
  print(original_plot)
  dev.off()
  cat("Original data plot saved to basicmotions_original_data.pdf\n")
  
  # Prepare data for WICMAD
  cat("\n3. Preparing data for WICMAD...\n")
  wicmad_data <- prepare_wicmad_data(data$train_series, data$train_labels)
  
  cat("Number of series:", length(wicmad_data$Y), "\n")
  cat("Time points:", length(wicmad_data$t), "\n")
  cat("Series dimensions:", nrow(wicmad_data$Y[[1]]), "x", ncol(wicmad_data$Y[[1]]), "\n")
  
  # Run WICMAD clustering
  cat("\n4. Running WICMAD clustering...\n")
  wicmad_result <- wicmad(
    Y = wicmad_data$Y,
    t = wicmad_data$t,
    n_iter = 10000,
    burn = 3000,
    thin = 1,
    warmup_iters = 500
  )
  
  # Extract cluster assignments
  cluster_assignments <- wicmad_result$Z[nrow(wicmad_result$Z), ]
  
  # Map WICMAD clusters to binary classification (normal vs anomaly)
  # Find the cluster that contains the most normal samples (label 0)
  normal_samples <- which(wicmad_data$labels == 0)
  cluster_counts <- table(cluster_assignments[normal_samples])
  normal_cluster <- as.numeric(names(cluster_counts)[which.max(cluster_counts)])
  
  # Map to binary: normal cluster = 0, all others = 1 (anomaly)
  binary_cluster_assignments <- ifelse(cluster_assignments == normal_cluster, 0, 1)
  
  cat("WICMAD found", length(unique(cluster_assignments)), "clusters\n")
  cat("Normal cluster:", normal_cluster, "\n")
  cat("Binary mapping: cluster", normal_cluster, "-> normal (0), others -> anomaly (1)\n")
  
  # Calculate metrics
  cat("\n5. Calculating clustering metrics...\n")
  metrics <- calculate_clustering_metrics(wicmad_data$labels, binary_cluster_assignments)
  
  cat("Clustering Performance:\n")
  cat("Macro Precision:", round(metrics$Macro_Precision, 4), "\n")
  cat("Macro Recall:", round(metrics$Macro_Recall, 4), "\n")
  cat("Macro Specificity:", round(metrics$Macro_Specificity, 4), "\n")
  cat("Macro F1 Score:", round(metrics$Macro_F1, 4), "\n")
  cat("Overall Accuracy:", round(metrics$Overall_Accuracy, 4), "\n")
  cat("Adjusted Rand Index:", round(metrics$ARI, 4), "\n")
  
  # Plot clustering results
  cat("\n6. Creating clustering results visualization...\n")
  clustering_plot <- plot_clustering_results(
    data$train_series, 
    wicmad_data$labels, 
    binary_cluster_assignments,
    "BasicMotions - Clustering Results"
  )
  
  # Save clustering results plot
  pdf("basicmotions_clustering_results.pdf", width = 12, height = 10)
  print(clustering_plot)
  dev.off()
  cat("Clustering results plot saved to basicmotions_clustering_results.pdf\n")
  
  # Print final summary
  cat("\n=== Analysis Complete ===\n")
  cat("Generated files:\n")
  cat("- basicmotions_original_data.pdf: Original multivariate time series (overlapped)\n")
  cat("- basicmotions_clustering_results.pdf: Clustered multivariate time series\n")
  
  cat("\nFinal Performance:\n")
  cat("Macro Precision:", round(metrics$Macro_Precision, 4), "\n")
  cat("Macro Recall:", round(metrics$Macro_Recall, 4), "\n")
  cat("Macro F1 Score:", round(metrics$Macro_F1, 4), "\n")
  cat("Overall Accuracy:", round(metrics$Overall_Accuracy, 4), "\n")
  cat("Adjusted Rand Index:", round(metrics$ARI, 4), "\n")
}

# Run the analysis
main()
