#!/usr/bin/env Rscript

# ArrowHead Dataset Analysis with WICMAD
# This script loads the ArrowHead dataset, runs WICMAD clustering, and visualizes results

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(cluster)
# library(caret)
library(mclust)  # For adjustedRandIndex
library(numDeriv)  # For computing derivatives

# Load WICMAD package
library(devtools)
devtools::load_all()

# Function to load ArrowHead dataset
load_arrowhead_data <- function(data_dir) {
  cat("Loading ArrowHead dataset...\n")
  
  # Load training data
  train_file <- file.path(data_dir, "ArrowHead_TRAIN.ts")
  test_file <- file.path(data_dir, "ArrowHead_TEST.ts")
  
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
  labels <- integer(n_series)
  
  for (i in seq_len(n_series)) {
    line <- data_lines[i]
    
    # Split by colon to separate data from label
    parts <- strsplit(line, ":")[[1]]
    if (length(parts) != 2) {
      stop("Invalid data format in line ", i, " - expected 2 parts, got ", length(parts))
    }
    
    # Parse the data part (comma-separated values)
    data_values <- as.numeric(strsplit(parts[1], ",")[[1]])
    
    # Parse the label part (integer)
    labels[i] <- as.integer(parts[2])
    
    # Store as a single column matrix (univariate time series)
    series_list[[i]] <- matrix(data_values, ncol = 1)
  }
  
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
      
      # Parse the data part
      data_values <- as.numeric(strsplit(parts[1], ",")[[1]])
      
      # Parse the label part
      test_labels[i] <- as.integer(parts[2])
      
      # Store as a single column matrix
      test_series_list[[i]] <- matrix(data_values, ncol = 1)
    }
  }
  
  # Combine train and test data
  all_series <- c(series_list, test_series_list)
  all_labels <- c(labels, test_labels)
  
  cat("Combined dataset:\n")
  cat("Train series:", length(series_list), "\n")
  cat("Test series:", length(test_series_list), "\n")
  cat("Total series:", length(all_series), "\n")
  
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
  
  # Create imbalanced dataset with approximately 5% anomalies
  normal_indices <- which(anomaly_labels == 0)
  anomaly_indices <- which(anomaly_labels == 1)
  
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
  cat("Target anomaly percentage: 5.0%\n")
  cat("Available normal samples used:", n_normal_used, "\n")
  cat("Available anomaly samples used:", n_anomalies_needed, "out of", length(anomaly_indices), "\n")
  
  return(list(
    train_series = imbalanced_series,
    train_labels = imbalanced_labels,
    test_series = NULL,
    test_labels = NULL
  ))
}

# Function to interpolate time series to target dimension
interpolate_series <- function(series, target_dim = 16) {
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

# Function to compute derivatives using numDeriv
compute_derivatives <- function(series) {
  n <- length(series)
  
  # Create interpolation function for smooth derivatives
  x_vals <- seq(0, 1, length.out = n)
  interp_func <- splinefun(x_vals, series, method = "natural")
  
  # Compute first derivative using numDeriv
  first_deriv <- sapply(x_vals, function(x) {
    tryCatch({
      grad(interp_func, x)
    }, error = function(e) {
      # Fallback to finite differences if grad fails
      if (x == 0) return((interp_func(x_vals[2]) - interp_func(x_vals[1])) / (x_vals[2] - x_vals[1]))
      if (x == 1) return((interp_func(x_vals[n]) - interp_func(x_vals[n-1])) / (x_vals[n] - x_vals[n-1]))
      h <- (x_vals[2] - x_vals[1])
      return((interp_func(x + h) - interp_func(x - h)) / (2 * h))
    })
  })
  
  # Compute second derivative using numDeriv
  second_deriv <- sapply(x_vals, function(x) {
    tryCatch({
      hessian(interp_func, x)
    }, error = function(e) {
      # Fallback to finite differences if hessian fails
      if (x == 0 || x == 1) return(0)
      h <- (x_vals[2] - x_vals[1])
      return((interp_func(x + h) - 2*interp_func(x) + interp_func(x - h)) / (h^2))
    })
  })
  
  return(list(first = first_deriv, second = second_deriv))
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
  # P=p_dim (time points), M=1 (univariate)
  Y <- interpolated_series
  
  return(list(Y = Y, t = t, labels = labels))
}

# Function to prepare derivative data for WICMAD
prepare_derivative_data <- function(series_list, labels) {
  cat("Preparing derivative data for WICMAD...\n")
  
  n_series <- length(series_list)
  
  # Find the highest power of 2 supported by the data
  max_dim <- max(sapply(series_list, nrow))
  p_dim <- 2^floor(log2(max_dim))
  cat("Using p parameter for derivatives (highest power of 2):", p_dim, "\n")
  
  # Interpolate original data to p dimensions
  original_interp <- matrix(0, nrow = n_series, ncol = p_dim)
  for (i in seq_len(n_series)) {
    original_interp[i, ] <- interpolate_series(series_list[[i]], p_dim)[, 1]  # Extract univariate data
  }
  
  # Compute derivatives for each interpolated series
  first_deriv_data <- matrix(0, nrow = n_series, ncol = p_dim)
  second_deriv_data <- matrix(0, nrow = n_series, ncol = p_dim)
  
  for (i in seq_len(n_series)) {
    derivs <- compute_derivatives(original_interp[i, ])
    first_deriv_data[i, ] <- derivs$first
    second_deriv_data[i, ] <- derivs$second
  }
  
  # Create time coordinates
  t <- seq(0, 1, length.out = p_dim)
  
  # Convert to list format for WICMAD (each element is a PxM matrix)
  # P=p_dim (time points), M=3 (original + first + second derivatives)
  Y <- lapply(seq_len(n_series), function(i) {
    matrix(c(original_interp[i, ], first_deriv_data[i, ], second_deriv_data[i, ]), nrow = p_dim, ncol = 3)
  })
  
  return(list(Y = Y, t = t, labels = labels, 
              original_data = original_interp,
              first_deriv_data = first_deriv_data,
              second_deriv_data = second_deriv_data))
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
plot_overlapped_data <- function(series_list, labels, title = "ArrowHead - Before Clustering", max_series = NULL) {
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
    
    series_df <- data.frame(
      time = time_points,
      value = series_data[, 1],  # Univariate data
      series_id = i,
      label = labels[i],
      label_name = ifelse(labels[i] == 0, "Normal", "Anomaly")
    )
    
    plot_data <- rbind(plot_data, series_df)
  }
  
  # Create overlapped plot
  p <- ggplot(plot_data, aes(x = .data$time, y = .data$value, color = .data$label_name, group = .data$series_id)) +
    geom_line(alpha = 0.6, linewidth = 0.5) +
    labs(title = title,
         x = "Time", y = "Value", color = "Class") +
    scale_color_manual(values = c("Normal" = "blue", "Anomaly" = "red")) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(p)
}

# Function to plot clustering results
plot_clustering_results <- function(series_list, true_labels, cluster_assignments, title = "ArrowHead - After Clustering", max_series = NULL) {
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
    
    series_df <- data.frame(
      time = time_points,
      value = series_data[, 1],  # Univariate data
      series_id = i,
      true_label = true_labels[i],
      cluster = cluster_assignments[i],
      cluster_name = ifelse(cluster_assignments[i] == 0, "Normal", "Anomaly")
    )
    
    plot_data <- rbind(plot_data, series_df)
  }
  
  # Create overlapped plot
  p <- ggplot(plot_data, aes(x = .data$time, y = .data$value, color = .data$cluster_name, group = .data$series_id)) +
    geom_line(alpha = 0.6, linewidth = 0.5) +
    labs(title = title,
         x = "Time", y = "Value", color = "Class") +
    scale_color_manual(values = c("Normal" = "blue", "Anomaly" = "red")) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(p)
}

# Function to plot derivative clustering results (3 separate channels)
plot_derivative_clustering_results <- function(original_data, first_deriv_data, second_deriv_data, cluster_assignments, title = "ArrowHead - After Clustering (Derivatives)", max_series = NULL) {
  cat("Creating derivative clustering results plot...\n")
  
  # Use all series if max_series is NULL
  if (is.null(max_series)) {
    series_indices <- seq_len(nrow(original_data))
  } else {
    n_series <- min(max_series, nrow(original_data))
    series_indices <- sample(seq_len(nrow(original_data)), n_series)
  }
  
  # Create plotting data for all three channels
  plot_data <- data.frame()
  
  for (i in series_indices) {
    time_points <- seq(0, 1, length.out = ncol(original_data))
    
    # Original signal
    orig_df <- data.frame(
      time = time_points,
      value = original_data[i, ],
      channel = "Original",
      series_id = i,
      cluster = cluster_assignments[i],
      cluster_name = ifelse(cluster_assignments[i] == 0, "Normal", "Anomaly")
    )
    
    # First derivative
    first_df <- data.frame(
      time = time_points,
      value = first_deriv_data[i, ],
      channel = "1st Derivative",
      series_id = i,
      cluster = cluster_assignments[i],
      cluster_name = ifelse(cluster_assignments[i] == 0, "Normal", "Anomaly")
    )
    
    # Second derivative
    second_df <- data.frame(
      time = time_points,
      value = second_deriv_data[i, ],
      channel = "2nd Derivative",
      series_id = i,
      cluster = cluster_assignments[i],
      cluster_name = ifelse(cluster_assignments[i] == 0, "Normal", "Anomaly")
    )
    
    plot_data <- rbind(plot_data, orig_df, first_df, second_df)
  }
  
  # Create faceted plot
  p <- ggplot(plot_data, aes(x = .data$time, y = .data$value, color = .data$cluster_name, group = .data$series_id)) +
    geom_line(alpha = 0.6, linewidth = 0.5) +
    facet_wrap(~ .data$channel, scales = "free_y", ncol = 1) +
    labs(title = title,
         x = "Time", y = "Value", color = "Class") +
    scale_color_manual(values = c("Normal" = "blue", "Anomaly" = "red")) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(p)
}

# Main analysis
main <- function() {
  cat("=== ArrowHead Dataset Analysis with WICMAD ===\n\n")
  
  # Ensure output directory exists
  output_dir <- "./plots/arrowhead"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("Created output directory:", output_dir, "\n")
  }
  
  # Set data directory
  data_dir <- "./data/ArrowHead"
  
  # Load data
  cat("1. Loading ArrowHead dataset...\n")
  data <- load_arrowhead_data(data_dir)
  
  # Use all observations in the normal class, then make final dataset 5% anomalies
  # (This is already handled in the load_arrowhead_data function)
  
  # Print data summary
  cat("Training data: ", length(data$train_series), " series\n")
  cat("Label distribution:\n")
  print(table(data$train_labels))
  cat("First few labels:", head(data$train_labels, 10), "\n")
  
  # Plot original data (overlapped)
  cat("\n2. Creating original data visualization...\n")
  cat("Number of series to plot:", length(data$train_series), "\n")
  cat("Number of labels:", length(data$train_labels), "\n")
  
  original_plot <- plot_overlapped_data(data$train_series, data$train_labels, 
                                       "ArrowHead Dataset - Before Clustering")
  cat("Plot object created successfully\n")
  
  # Save original data plot
  plot_path <- "./plots/arrowhead/arrowhead_original_data.pdf"
  cat("Attempting to save plot to:", plot_path, "\n")
  cat("Current working directory:", getwd(), "\n")
  
  pdf(plot_path, width = 12, height = 8)
  print(original_plot)
  dev.off()
  
  # Check if file was created
  if (file.exists(plot_path)) {
    cat("✓ Original data plot saved successfully to:", plot_path, "\n")
  } else {
    cat("✗ ERROR: Plot file was not created at:", plot_path, "\n")
  }
  
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
    n_iter = 8000,
    burn = 3000,
    warmup_iters = 500,
    unpin = FALSE
  )
  
  # Extract cluster assignments
  cluster_assignments <- wicmad_result$Z[nrow(wicmad_result$Z), ]
  
  cat("Raw cluster assignments:", cluster_assignments, "\n")
  cat("Unique clusters found:", unique(cluster_assignments), "\n")
  cat("True labels:", wicmad_data$labels, "\n")
  
  # Map WICMAD clusters to binary classification (normal vs anomaly)
  # Find the cluster that contains the most normal samples (label 0)
  normal_samples <- which(wicmad_data$labels == 0)
  cat("Normal sample indices:", normal_samples, "\n")
  cat("Normal sample cluster assignments:", cluster_assignments[normal_samples], "\n")
  
  cluster_counts <- table(cluster_assignments[normal_samples])
  cat("Cluster counts for normal samples:", cluster_counts, "\n")
  
  if (length(cluster_counts) > 0) {
    normal_cluster <- as.numeric(names(cluster_counts)[which.max(cluster_counts)])
  } else {
    # If no normal samples, use the most common cluster overall
    overall_counts <- table(cluster_assignments)
    normal_cluster <- as.numeric(names(overall_counts)[which.max(overall_counts)])
  }
  
  # Map to binary: normal cluster = 0, all others = 1 (anomaly)
  binary_cluster_assignments <- ifelse(cluster_assignments == normal_cluster, 0, 1)
  
  cat("WICMAD found", length(unique(cluster_assignments)), "clusters\n")
  cat("Normal cluster:", normal_cluster, "\n")
  cat("Binary mapping: cluster", normal_cluster, "-> normal (0), others -> anomaly (1)\n")
  cat("Binary cluster assignments:", binary_cluster_assignments, "\n")
  
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
    "ArrowHead - After Clustering"
  )
  
  # Save clustering results plot
  clustering_plot_path <- "./plots/arrowhead/arrowhead_clustering_results.pdf"
  cat("Attempting to save clustering plot to:", clustering_plot_path, "\n")
  
  pdf(clustering_plot_path, width = 12, height = 8)
  print(clustering_plot)
  dev.off()
  
  # Check if file was created
  if (file.exists(clustering_plot_path)) {
    cat("✓ Clustering results plot saved successfully to:", clustering_plot_path, "\n")
  } else {
    cat("✗ ERROR: Clustering plot file was not created at:", clustering_plot_path, "\n")
  }
  
  # Run analysis on derivatives
  cat("\n7. Analyzing derivative data...\n")
  deriv_wicmad_data <- prepare_derivative_data(data$train_series, data$train_labels)
  
  cat("Number of series:", length(deriv_wicmad_data$Y), "\n")
  cat("Time points:", length(deriv_wicmad_data$t), "\n")
  cat("Series dimensions:", nrow(deriv_wicmad_data$Y[[1]]), "x", ncol(deriv_wicmad_data$Y[[1]]), "\n")
  
  # Run WICMAD clustering on derivatives
  cat("Running WICMAD clustering on derivatives...\n")
  deriv_wicmad_result <- wicmad(
    Y = deriv_wicmad_data$Y,
    t = deriv_wicmad_data$t,
    n_iter = 8000,
    burn = 3000,
    warmup_iters = 500,
    unpin = FALSE
  )
  
  # Extract cluster assignments
  deriv_cluster_assignments <- deriv_wicmad_result$Z[nrow(deriv_wicmad_result$Z), ]
  
  cat("Raw derivative cluster assignments:", deriv_cluster_assignments, "\n")
  cat("Unique derivative clusters found:", unique(deriv_cluster_assignments), "\n")
  
  # Map WICMAD clusters to binary classification (normal vs anomaly)
  normal_samples <- which(deriv_wicmad_data$labels == 0)
  deriv_cluster_counts <- table(deriv_cluster_assignments[normal_samples])
  cat("Derivative cluster counts for normal samples:", deriv_cluster_counts, "\n")
  
  if (length(deriv_cluster_counts) > 0) {
    deriv_normal_cluster <- as.numeric(names(deriv_cluster_counts)[which.max(deriv_cluster_counts)])
  } else {
    overall_counts <- table(deriv_cluster_assignments)
    deriv_normal_cluster <- as.numeric(names(overall_counts)[which.max(overall_counts)])
  }
  
  # Map to binary: normal cluster = 0, all others = 1 (anomaly)
  deriv_binary_cluster_assignments <- ifelse(deriv_cluster_assignments == deriv_normal_cluster, 0, 1)
  
  cat("WICMAD found", length(unique(deriv_cluster_assignments)), "derivative clusters\n")
  cat("Derivative normal cluster:", deriv_normal_cluster, "\n")
  cat("Derivative binary cluster assignments:", deriv_binary_cluster_assignments, "\n")
  
  # Calculate metrics for derivatives
  cat("Calculating derivative clustering metrics...\n")
  deriv_metrics <- calculate_clustering_metrics(deriv_wicmad_data$labels, deriv_binary_cluster_assignments)
  
  cat("Derivative Clustering Performance:\n")
  cat("Macro Precision:", round(deriv_metrics$Macro_Precision, 4), "\n")
  cat("Macro Recall:", round(deriv_metrics$Macro_Recall, 4), "\n")
  cat("Macro Specificity:", round(deriv_metrics$Macro_Specificity, 4), "\n")
  cat("Macro F1 Score:", round(deriv_metrics$Macro_F1, 4), "\n")
  cat("Overall Accuracy:", round(deriv_metrics$Overall_Accuracy, 4), "\n")
  cat("Adjusted Rand Index:", round(deriv_metrics$ARI, 4), "\n")
  
  # Plot derivative clustering results (3 separate channels)
  deriv_clustering_plot <- plot_derivative_clustering_results(
    deriv_wicmad_data$original_data,
    deriv_wicmad_data$first_deriv_data,
    deriv_wicmad_data$second_deriv_data,
    deriv_binary_cluster_assignments,
    "ArrowHead - After Clustering (Derivatives)"
  )
  
  # Save derivative clustering results plot
  deriv_clustering_plot_path <- "./plots/arrowhead/arrowhead_derivatives_clustering.pdf"
  cat("Attempting to save derivative clustering plot to:", deriv_clustering_plot_path, "\n")
  
  pdf(deriv_clustering_plot_path, width = 12, height = 8)
  print(deriv_clustering_plot)
  dev.off()
  
  # Check if file was created
  if (file.exists(deriv_clustering_plot_path)) {
    cat("✓ Derivative clustering results plot saved successfully to:", deriv_clustering_plot_path, "\n")
  } else {
    cat("✗ ERROR: Derivative clustering plot file was not created at:", deriv_clustering_plot_path, "\n")
  }
  
  # Print final summary
  cat("\n=== Analysis Complete ===\n")
  cat("Generated files:\n")
  cat("- ./plots/arrowhead/arrowhead_original_data.pdf: Original time series (overlapped)\n")
  cat("- ./plots/arrowhead/arrowhead_clustering_results.pdf: Clustered time series\n")
  cat("- ./plots/arrowhead/arrowhead_derivatives_clustering.pdf: Derivative clustering results\n")
  
  cat("\nPerformance Comparison:\n")
  cat("Raw Signal - Precision:", round(metrics$Macro_Precision, 4), 
      "Recall:", round(metrics$Macro_Recall, 4),
      "F1:", round(metrics$Macro_F1, 4), 
      "Accuracy:", round(metrics$Overall_Accuracy, 4),
      "ARI:", round(metrics$ARI, 4), "\n")
  cat("Derivatives - Precision:", round(deriv_metrics$Macro_Precision, 4), 
      "Recall:", round(deriv_metrics$Macro_Recall, 4),
      "F1:", round(deriv_metrics$Macro_F1, 4), 
      "Accuracy:", round(deriv_metrics$Overall_Accuracy, 4),
      "ARI:", round(deriv_metrics$ARI, 4), "\n")
}

# Run the analysis
main()
