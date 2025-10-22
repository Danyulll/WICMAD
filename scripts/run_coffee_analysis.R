#!/usr/bin/env Rscript

# Coffee Dataset Analysis with WICMAD
# This script loads the coffee dataset, runs WICMAD clustering, and visualizes results

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(cluster)
# library(caret)
library(fda)  # For FPCA
library(mclust)  # For adjustedRandIndex

# Load WICMAD package
library(devtools)
devtools::install()

# Function to load coffee dataset and create imbalanced data
load_coffee_data <- function(data_dir, imbalance_ratio = 0.05) {
  cat("Loading coffee dataset...\n")
  
  # Load training data
  train_file <- file.path(data_dir, "Coffee_TRAIN.txt")
  test_file <- file.path(data_dir, "Coffee_TEST.txt")
  
  if (!file.exists(train_file)) {
    stop("Training file not found: ", train_file)
  }
  
  # Load training data
  train_data <- read.table(train_file, header = FALSE)
  
  # Load test data if available
  test_data <- NULL
  if (file.exists(test_file)) {
    test_data <- read.table(test_file, header = FALSE)
  }
  
  # Extract labels (first column) and features (remaining columns)
  train_labels <- train_data[, 1]
  train_features <- as.matrix(train_data[, -1])
  
  test_labels <- NULL
  test_features <- NULL
  if (!is.null(test_data)) {
    test_labels <- test_data[, 1]
    test_features <- as.matrix(test_data[, -1])
  }
  
  # If no test data, create synthetic data to have enough samples for 5% imbalance
  if (is.null(test_data)) {
    cat("No test data found. Creating synthetic data for proper imbalance...\n")
    # Duplicate the training data to have more samples
    n_duplicates <- 3  # This will give us 28 * 4 = 112 total samples
    train_labels <- rep(train_labels, n_duplicates + 1)
    train_features <- rbind(train_features, 
                           train_features[rep(seq_len(nrow(train_features)), n_duplicates), ])
    test_labels <- NULL
    test_features <- NULL
  }
  
  # Create imbalanced dataset
  cat("Creating imbalanced dataset with", round(imbalance_ratio * 100, 1), "% anomaly class...\n")
  
  # Combine train and test data for resampling
  all_labels <- c(train_labels, test_labels)
  all_features <- rbind(train_features, test_features)
  
  # Find indices for each class
  normal_indices <- which(all_labels == 0)
  anomaly_indices <- which(all_labels == 1)
  
  cat("Original distribution:\n")
  cat("Normal samples:", length(normal_indices), "\n")
  cat("Anomaly samples:", length(anomaly_indices), "\n")
  
  # Calculate target numbers for imbalanced dataset
  total_samples <- length(all_labels)
  target_anomaly <- max(1, round(total_samples * imbalance_ratio))
  target_normal <- total_samples - target_anomaly
  
  cat("Target distribution:\n")
  cat("Normal samples:", target_normal, "\n")
  cat("Anomaly samples:", target_anomaly, "\n")
  
  # Ensure we have enough samples of each class
  if (length(anomaly_indices) < target_anomaly) {
    # If we don't have enough anomalies, duplicate some
    needed_anomaly <- target_anomaly - length(anomaly_indices)
    additional_anomaly <- sample(anomaly_indices, needed_anomaly, replace = TRUE)
    selected_anomaly <- c(anomaly_indices, additional_anomaly)
  } else {
    selected_anomaly <- sample(anomaly_indices, target_anomaly)
  }
  
  # Sample normal indices
  if (length(normal_indices) >= target_normal) {
    selected_normal <- sample(normal_indices, target_normal)
  } else {
    # If we don't have enough normal samples, duplicate some
    needed_normal <- target_normal - length(normal_indices)
    additional_normal <- sample(normal_indices, needed_normal, replace = TRUE)
    selected_normal <- c(normal_indices, additional_normal)
  }
  
  # Combine selected indices
  selected_indices <- c(selected_normal, selected_anomaly)
  
  # Create imbalanced dataset
  imbalanced_labels <- all_labels[selected_indices]
  imbalanced_features <- all_features[selected_indices, ]
  
  # Shuffle the data
  shuffle_indices <- sample(length(imbalanced_labels))
  imbalanced_labels <- imbalanced_labels[shuffle_indices]
  imbalanced_features <- imbalanced_features[shuffle_indices, ]
  
  cat("Final imbalanced distribution:\n")
  print(table(imbalanced_labels))
  cat("Anomaly percentage:", round(mean(imbalanced_labels == 1) * 100, 1), "%\n")
  
  return(list(
    train_labels = imbalanced_labels,
    train_features = imbalanced_features,
    test_labels = NULL,
    test_features = NULL
  ))
}

# Function to interpolate time series to target dimension
interpolate_series <- function(series, target_dim = 128) {
  n <- length(series)
  if (n == target_dim) {
    return(series)
  }
  
  # Create time points for interpolation
  x_old <- seq(0, 1, length.out = n)
  x_new <- seq(0, 1, length.out = target_dim)
  
  # Interpolate
  interpolated <- approx(x_old, series, x_new, method = "linear")$y
  
  return(interpolated)
}

# Function to compute derivatives
compute_derivatives <- function(series) {
  # First derivative
  first_deriv <- diff(series)
  
  # Second derivative
  second_deriv <- diff(first_deriv)
  
  # Pad to same length as original
  first_deriv <- c(first_deriv[1], first_deriv)
  second_deriv <- c(second_deriv[1], second_deriv[1], second_deriv)
  
  return(list(first = first_deriv, second = second_deriv))
}

# Function to perform functional PCA and select components for 90% variance
perform_fpca <- function(features, target_dim = 128) {
  cat("Performing functional PCA analysis...\n")
  
  # Interpolate all series to same length
  n_series <- nrow(features)
  interpolated_features <- matrix(0, nrow = n_series, ncol = target_dim)
  
  for (i in seq_len(n_series)) {
    interpolated_features[i, ] <- interpolate_series(features[i, ], target_dim)
  }
  
  # Create time points
  t <- seq(0, 1, length.out = target_dim)
  
  # Perform functional PCA using fda package
  # Create functional data object
  fd_obj <- fda::Data2fd(argvals = t, y = t(interpolated_features))
  
  # Perform functional PCA
  fpca_result <- fda::pca.fd(fd_obj, nharm = min(10, n_series-1))
  
  # Calculate cumulative variance explained
  cumvar <- cumsum(fpca_result$varprop)
  
  # Find number of components for 90% variance
  n_components <- which(cumvar >= 0.9)[1]
  
  # Ensure n_components is dyadic (power of 2) for WICMAD
  dyadic_components <- 2^ceiling(log2(n_components))
  if (dyadic_components > length(fpca_result$varprop)) {
    dyadic_components <- 2^floor(log2(length(fpca_result$varprop)))
  }
  
  cat("Number of components for 90% variance:", n_components, "\n")
  cat("Actual variance explained:", round(cumvar[n_components], 4), "\n")
  cat("Using dyadic components:", dyadic_components, "\n")
  cat("Dyadic variance explained:", round(cumvar[dyadic_components], 4), "\n")
  
  n_components <- dyadic_components
  
  # Get functional principal component scores
  fpca_scores <- fpca_result$scores[, 1:n_components, drop = FALSE]
  
  # Get functional principal component functions (eigenfunctions)
  fpca_functions <- fpca_result$harmonics[1:n_components]
  
  return(list(
    scores = fpca_scores,
    functions = fpca_functions,
    n_components = n_components,
    cumvar = cumvar[n_components],
    fpca_result = fpca_result,
    time_points = t
  ))
}

# Function to prepare data for WICMAD
prepare_wicmad_data <- function(features, labels) {
  cat("Preparing data for WICMAD...\n")
  
  # Interpolate each time series to 32 dimensions
  n_series <- nrow(features)
  interpolated_features <- matrix(0, nrow = n_series, ncol = 32)
  
  for (i in seq_len(n_series)) {
    interpolated_features[i, ] <- interpolate_series(features[i, ], 32)
  }
  
  # Create time coordinates (normalized to [0,1])
  t <- seq(0, 1, length.out = 32)
  
  # Convert to list format for WICMAD (each element is a PxM matrix)
  # For univariate data, we need P=32 (time points) and M=1 (channels)
  Y <- lapply(seq_len(n_series), function(i) {
    matrix(interpolated_features[i, ], nrow = 32, ncol = 1)
  })
  
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
plot_overlapped_data <- function(features, labels, title = "Coffee Dataset", max_series = NULL) {
  cat("Creating overlapped data plot...\n")
  
  # Use all series if max_series is NULL
  if (is.null(max_series)) {
    series_indices <- seq_len(nrow(features))
  } else {
    n_series <- min(max_series, nrow(features))
    series_indices <- sample(seq_len(nrow(features)), n_series)
  }
  
  # Create plotting data
  plot_data <- data.frame()
  
  for (i in series_indices) {
    series_data <- features[i, ]
    time_points <- seq(0, 1, length.out = length(series_data))
    
    series_df <- data.frame(
      time = time_points,
      value = series_data,
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

# Function to plot original data (overlapped)
plot_original_data <- function(features, labels, max_series = NULL) {
  return(plot_overlapped_data(features, labels, "Coffee Dataset - Original Time Series", max_series))
}

# Function to plot clustering results (overlapped)
plot_clustering_results <- function(features, true_labels, cluster_assignments, title = "Coffee Dataset - Clustering Results", max_series = NULL) {
  cat("Creating clustering results plot...\n")
  
  # Use all series if max_series is NULL
  if (is.null(max_series)) {
    series_indices <- seq_len(nrow(features))
  } else {
    n_series <- min(max_series, nrow(features))
    series_indices <- sample(seq_len(nrow(features)), n_series)
  }
  
  # Create plotting data
  plot_data <- data.frame()
  
  for (i in series_indices) {
    series_data <- features[i, ]
    time_points <- seq(0, 1, length.out = length(series_data))
    
    series_df <- data.frame(
      time = time_points,
      value = series_data,
      series_id = i,
      true_label = true_labels[i],
      cluster = cluster_assignments[i],
      true_label_name = ifelse(true_labels[i] == 0, "Normal", "Anomaly")
    )
    
    plot_data <- rbind(plot_data, series_df)
  }
  
  # Create overlapped plot
  p <- ggplot(plot_data, aes(x = .data$time, y = .data$value, color = .data$true_label_name, group = .data$series_id)) +
    geom_line(alpha = 0.6, linewidth = 0.5) +
    labs(title = title,
         x = "Time", y = "Value", color = "True Class") +
    scale_color_manual(values = c("Normal" = "blue", "Anomaly" = "red")) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(p)
}

# Function to plot derivative clustering results (3 separate channels)
plot_derivative_clustering_results <- function(original_data, first_deriv_data, second_deriv_data, true_labels, cluster_assignments, title = "Derivatives - Clustering Results", max_series = NULL) {
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
      true_label = true_labels[i],
      cluster = cluster_assignments[i],
      true_label_name = ifelse(true_labels[i] == 0, "Normal", "Anomaly")
    )
    
    # First derivative
    first_df <- data.frame(
      time = time_points,
      value = first_deriv_data[i, ],
      channel = "1st Derivative",
      series_id = i,
      true_label = true_labels[i],
      cluster = cluster_assignments[i],
      true_label_name = ifelse(true_labels[i] == 0, "Normal", "Anomaly")
    )
    
    # Second derivative
    second_df <- data.frame(
      time = time_points,
      value = second_deriv_data[i, ],
      channel = "2nd Derivative",
      series_id = i,
      true_label = true_labels[i],
      cluster = cluster_assignments[i],
      true_label_name = ifelse(true_labels[i] == 0, "Normal", "Anomaly")
    )
    
    plot_data <- rbind(plot_data, orig_df, first_df, second_df)
  }
  
  # Create faceted plot
  p <- ggplot(plot_data, aes(x = .data$time, y = .data$value, color = .data$true_label_name, group = .data$series_id)) +
    geom_line(alpha = 0.6, linewidth = 0.5) +
    facet_wrap(~ .data$channel, scales = "free_y", ncol = 1) +
    labs(title = title,
         x = "Time", y = "Value", color = "True Class") +
    scale_color_manual(values = c("Normal" = "blue", "Anomaly" = "red")) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(p)
}

# Function to plot FPCA clustering results (showing functional components)
plot_fpca_clustering_results <- function(fpca_result, true_labels, cluster_assignments, title = "FPCA - Clustering Results", max_series = NULL) {
  cat("Creating FPCA clustering results plot...\n")
  
  # Use all series if max_series is NULL
  n_series <- nrow(fpca_result$scores)
  if (is.null(max_series)) {
    series_indices <- seq_len(n_series)
  } else {
    n_series <- min(max_series, n_series)
    series_indices <- sample(seq_len(n_series), n_series)
  }
  
  # Create plotting data for functional components
  plot_data <- data.frame()
  
  for (i in series_indices) {
    time_points <- fpca_result$time_points
    
    for (j in seq_len(fpca_result$n_components)) {
      # Evaluate the j-th functional component at time points
      component_values <- fda::eval.fd(time_points, fpca_result$functions[j])
      
      component_df <- data.frame(
        time = time_points,
        value = as.vector(component_values),
        component = j,
        component_name = paste("Component", j),
        series_id = i,
        true_label = true_labels[i],
        cluster = cluster_assignments[i],
        true_label_name = ifelse(true_labels[i] == 0, "Normal", "Anomaly"),
        stringsAsFactors = FALSE
      )
      
      plot_data <- rbind(plot_data, component_df)
    }
  }
  
  # Create faceted plot for functional components
  p <- ggplot(plot_data, aes(x = .data$time, y = .data$value, color = .data$true_label_name, group = .data$series_id)) +
    geom_line(alpha = 0.6, linewidth = 0.5) +
    facet_wrap(~ .data$component_name, scales = "free_y", ncol = 2) +
    labs(title = title,
         x = "Time", y = "Component Value", color = "True Class") +
    scale_color_manual(values = c("Normal" = "blue", "Anomaly" = "red")) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(p)
}

# Function to prepare derivative data for WICMAD
prepare_derivative_data <- function(features, labels) {
  cat("Preparing derivative data for WICMAD...\n")
  
  n_series <- nrow(features)
  
  # Interpolate original data to 128 dimensions first
  original_interp <- matrix(0, nrow = n_series, ncol = 128)
  for (i in seq_len(n_series)) {
    original_interp[i, ] <- interpolate_series(features[i, ], 128)
  }
  
  # Compute derivatives for each interpolated series
  first_deriv_data <- matrix(0, nrow = n_series, ncol = 128)
  second_deriv_data <- matrix(0, nrow = n_series, ncol = 128)
  
  for (i in seq_len(n_series)) {
    derivs <- compute_derivatives(original_interp[i, ])
    first_deriv_data[i, ] <- derivs$first
    second_deriv_data[i, ] <- derivs$second
  }
  
  # Create time coordinates
  t <- seq(0, 1, length.out = 128)
  
  # Convert to list format for WICMAD (each element is a PxM matrix)
  # P=128 (time points), M=3 (original + first + second derivatives)
  Y <- lapply(seq_len(n_series), function(i) {
    matrix(c(original_interp[i, ], first_deriv_data[i, ], second_deriv_data[i, ]), nrow = 128, ncol = 3)
  })
  
  return(list(Y = Y, t = t, labels = labels, 
              original_data = original_interp,
              first_deriv_data = first_deriv_data,
              second_deriv_data = second_deriv_data))
}

# Function to prepare FPCA data for WICMAD
prepare_fpca_data <- function(features, labels) {
  cat("Preparing FPCA data for WICMAD...\n")
  
  # Perform functional PCA
  fpca_result <- perform_fpca(features, 128)
  
  # Create time coordinates for the functional components
  t <- fpca_result$time_points
  
  # Convert to list format for WICMAD
  # Each element is a PxM matrix where P=time_points, M=n_components (functional components)
  Y <- lapply(seq_len(nrow(fpca_result$scores)), function(i) {
    # Get the functional components for this series
    fpca_scores <- fpca_result$scores[i, ]
    
    # Create a matrix where each column is a functional component
    # We need to evaluate the functional components at the time points
    component_matrix <- matrix(0, nrow = length(t), ncol = fpca_result$n_components)
    
    for (j in seq_len(fpca_result$n_components)) {
      # Evaluate the j-th functional component at time points
      component_values <- fda::eval.fd(t, fpca_result$functions[j])
      component_matrix[, j] <- component_values
    }
    
    component_matrix
  })
  
  return(list(Y = Y, t = t, labels = labels, fpca_result = fpca_result))
}

# Function to run WICMAD analysis on a dataset
run_wicmad_analysis <- function(data_prep, data_name, features, labels) {
  cat(paste("\n=== Running WICMAD on", data_name, "===\n"))
  
  # Prepare data
  wicmad_data <- data_prep(features, labels)
  
  cat("Number of series:", length(wicmad_data$Y), "\n")
  cat("Time points:", length(wicmad_data$t), "\n")
  cat("Series dimensions:", nrow(wicmad_data$Y[[1]]), "x", ncol(wicmad_data$Y[[1]]), "\n")
  
  # Run WICMAD
  cat("Running WICMAD clustering...\n")
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
  
  # Calculate metrics
  metrics <- calculate_clustering_metrics(wicmad_data$labels, cluster_assignments)
  
  cat("Clustering Performance:\n")
  cat("Macro Precision:", round(metrics$Macro_Precision, 4), "\n")
  cat("Macro Recall:", round(metrics$Macro_Recall, 4), "\n")
  cat("Macro Specificity:", round(metrics$Macro_Specificity, 4), "\n")
  cat("Macro F1 Score:", round(metrics$Macro_F1, 4), "\n")
  cat("Overall Accuracy:", round(metrics$Overall_Accuracy, 4), "\n")
  cat("Adjusted Rand Index:", round(metrics$ARI, 4), "\n")
  
  return(list(
    metrics = metrics,
    cluster_assignments = cluster_assignments,
    wicmad_result = wicmad_result
  ))
}

# Main analysis
main <- function() {
  cat("=== Coffee Dataset Analysis with WICMAD ===\n\n")
  
  # Set data directory
  data_dir <- "data/Coffee"
  
  # Load data with 5% anomaly class
  cat("1. Loading coffee dataset with 5% anomaly class...\n")
  data <- load_coffee_data(data_dir, imbalance_ratio = 0.05)
  
  # Print data summary
  cat("Training data shape:", nrow(data$train_features), "x", ncol(data$train_features), "\n")
  cat("Label distribution:\n")
  print(table(data$train_labels))
  cat("First few labels:", head(data$train_labels, 10), "\n")
  
  # Plot original data (overlapped)
  cat("\n2. Creating original data visualization...\n")
  original_plot <- plot_overlapped_data(data$train_features, data$train_labels, 
                                       "Coffee Dataset - Imbalanced (5% Anomaly)")
  
  # Save original data plot
  pdf("plots/coffee_imbalanced_original_data.pdf", width = 12, height = 8)
  print(original_plot)
  dev.off()
  cat("Imbalanced original data plot saved to plots/coffee_imbalanced_original_data.pdf\n")
  
  # Run analysis on raw signal
  cat("\n3. Analyzing raw signal data...\n")
  raw_results <- run_wicmad_analysis(
    prepare_wicmad_data, 
    "Raw Signal", 
    data$train_features, 
    data$train_labels
  )
  
  # Plot raw signal clustering results
  raw_clustering_plot <- plot_clustering_results(
    data$train_features, 
    data$train_labels, 
    raw_results$cluster_assignments,
    "Raw Signal - Clustering Results"
  )
  
  pdf("plots/coffee_raw_signal_clustering.pdf", width = 12, height = 8)
  print(raw_clustering_plot)
  dev.off()
  cat("Raw signal clustering plot saved to plots/coffee_raw_signal_clustering.pdf\n")
  
  # Run analysis on derivatives
  cat("\n4. Analyzing derivative data...\n")
  deriv_results <- run_wicmad_analysis(
    prepare_derivative_data, 
    "Derivatives", 
    data$train_features, 
    data$train_labels
  )
  
  # Plot derivative clustering results (3 separate channels)
  deriv_data_prep <- prepare_derivative_data(data$train_features, data$train_labels)
  deriv_clustering_plot <- plot_derivative_clustering_results(
    deriv_data_prep$original_data,
    deriv_data_prep$first_deriv_data,
    deriv_data_prep$second_deriv_data,
    data$train_labels, 
    deriv_results$cluster_assignments,
    "Derivatives - Clustering Results"
  )
  
  pdf("plots/coffee_derivatives_clustering.pdf", width = 12, height = 8)
  print(deriv_clustering_plot)
  dev.off()
  cat("Derivatives clustering plot saved to plots/coffee_derivatives_clustering.pdf\n")
  
  # Run analysis on FPCA data
  cat("\n5. Analyzing FPCA data...\n")
  fpca_results <- run_wicmad_analysis(
    prepare_fpca_data, 
    "FPCA", 
    data$train_features, 
    data$train_labels
  )
  
  # Plot FPCA clustering results (showing functional components)
  # Get FPCA data from the preparation step
  fpca_data_prep <- prepare_fpca_data(data$train_features, data$train_labels)
  fpca_clustering_plot <- plot_fpca_clustering_results(
    fpca_data_prep$fpca_result, 
    data$train_labels, 
    fpca_results$cluster_assignments,
    "FPCA - Clustering Results"
  )
  
  pdf("plots/coffee_fpca_clustering.pdf", width = 12, height = 8)
  print(fpca_clustering_plot)
  dev.off()
  cat("FPCA clustering plot saved to plots/coffee_fpca_clustering.pdf\n")
  
  # Print final summary
  cat("\n=== Imbalanced Analysis Complete ===\n")
  cat("Generated files:\n")
  cat("- plots/coffee_imbalanced_original_data.pdf: Imbalanced time series (5% anomaly, overlapped)\n")
  cat("- plots/coffee_raw_signal_clustering.pdf: Raw signal clustering results\n")
  cat("- plots/coffee_derivatives_clustering.pdf: Derivatives clustering results\n")
  cat("- plots/coffee_fpca_clustering.pdf: FPCA clustering results\n")
  
  cat("\nPerformance Comparison:\n")
  cat("Raw Signal - Precision:", round(raw_results$metrics$Macro_Precision, 4), 
      "Recall:", round(raw_results$metrics$Macro_Recall, 4),
      "F1:", round(raw_results$metrics$Macro_F1, 4), 
      "Accuracy:", round(raw_results$metrics$Overall_Accuracy, 4),
      "ARI:", round(raw_results$metrics$ARI, 4), "\n")
  cat("Derivatives - Precision:", round(deriv_results$metrics$Macro_Precision, 4), 
      "Recall:", round(deriv_results$metrics$Macro_Recall, 4),
      "F1:", round(deriv_results$metrics$Macro_F1, 4), 
      "Accuracy:", round(deriv_results$metrics$Overall_Accuracy, 4),
      "ARI:", round(deriv_results$metrics$ARI, 4), "\n")
  cat("FPCA - Precision:", round(fpca_results$metrics$Macro_Precision, 4), 
      "Recall:", round(fpca_results$metrics$Macro_Recall, 4),
      "F1:", round(fpca_results$metrics$Macro_F1, 4), 
      "Accuracy:", round(fpca_results$metrics$Overall_Accuracy, 4),
      "ARI:", round(fpca_results$metrics$ARI, 4), "\n")
}

main()
