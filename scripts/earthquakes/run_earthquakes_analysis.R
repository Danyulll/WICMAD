#!/usr/bin/env Rscript

# Set working directory to scripts folder
# setwd("scripts")

# Load required libraries
library(devtools)
library(ggplot2)
library(dplyr)
library(WICMAD)
library(numDeriv)

# Load the WICMAD package
devtools::load_all()

# Function to load earthquakes dataset and create imbalanced data
load_earthquakes_data <- function(data_dir, imbalance_ratio = 0.05) {
  cat("Loading earthquakes dataset...\n")
  
  # Read training and test data
  train_file <- file.path(data_dir, "Earthquakes_TRAIN.txt")
  test_file <- file.path(data_dir, "Earthquakes_TEST.txt")
  
  # Load training data
  train_data <- read.table(train_file, header = FALSE)
  train_labels <- train_data[, 1]  # First column is labels
  train_series <- train_data[, -1]  # Rest are time series data
  
  # Load test data
  test_data <- read.table(test_file, header = FALSE)
  test_labels <- test_data[, 1]
  test_series <- test_data[, -1]
  
  # Combine train and test data
  all_series <- rbind(train_series, test_series)
  all_labels <- c(train_labels, test_labels)
  
  cat("Combined dataset size:", nrow(all_series), "x", ncol(all_series), "\n")
  cat("Original label distribution:\n")
  print(table(all_labels))
  
  # Determine majority class (0 = normal, 1 = earthquake)
  label_counts <- table(all_labels)
  majority_class <- as.numeric(names(label_counts)[which.max(label_counts)])
  minority_class <- 1 - majority_class
  
  cat("Majority class (normal):", majority_class, "\n")
  cat("Minority class (earthquake):", minority_class, "\n")
  
  # Create imbalanced dataset with 5% minority class
  # Use all available data for analysis
  cat("Using all available data for analysis\n")
  
  # Use the data as-is
  imbalanced_series <- all_series
  imbalanced_labels <- all_labels
  
  # Map to binary classification (0 = normal, 1 = anomaly)
  binary_labels <- ifelse(imbalanced_labels == majority_class, 0, 1)
  
  cat("Final dataset size:", nrow(imbalanced_series), "\n")
  cat("Final label distribution:\n")
  print(table(binary_labels))
  
  # Find the highest power of 2 supported by the data
  max_dim <- ncol(all_series)
  target_dim <- 2^floor(log2(max_dim))
  cat("Using p parameter (highest power of 2):", target_dim, "\n")
  
  return(list(
    train_data = imbalanced_series,
    train_labels = binary_labels,
    majority_class = 0,  # Normal is the majority class
    target_dim = target_dim
  ))
}

# Function to interpolate time series to target dimension
interpolate_series <- function(series, target_dim = 16) {
  n <- length(series)
  if (n == target_dim) {
    return(series)
  }
  
  # Create interpolation points
  x_old <- seq(0, 1, length.out = n)
  x_new <- seq(0, 1, length.out = target_dim)
  
  # Interpolate
  interpolated <- approx(x_old, series, x_new, method = "linear")$y
  
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


# Function to prepare data for WICMAD (raw time series)
prepare_wicmad_data <- function(series_data, labels, target_dim = 16) {
  cat("Preparing raw time series data for WICMAD...\n")
  
  # Interpolate each series to target dimensions
  n_series <- nrow(series_data)
  interpolated_series <- matrix(0, nrow = n_series, ncol = target_dim)
  
  for (i in seq_len(n_series)) {
    interpolated_series[i, ] <- interpolate_series(series_data[i, ], target_dim)
  }
  
  # Create time coordinates
  t <- seq(0, 1, length.out = target_dim)
  
  # Convert to list format for WICMAD (each element is a PxM matrix)
  # For univariate data, we need P=target_dim (time points) and M=1 (channels)
  Y <- lapply(seq_len(n_series), function(i) {
    matrix(interpolated_series[i, ], nrow = target_dim, ncol = 1)
  })
  
  return(list(Y = Y, t = t, labels = labels))
}

# Function to prepare derivative data for WICMAD
prepare_derivative_data <- function(series_data, labels, target_dim = 16) {
  cat("Preparing derivative data for WICMAD...\n")
  
  n_series <- nrow(series_data)
  
  # Compute derivatives for each series
  first_deriv_list <- list()
  second_deriv_list <- list()
  
  for (i in seq_len(n_series)) {
    # Interpolate to target dimension first
    interpolated <- interpolate_series(series_data[i, ], target_dim)
    
    # Compute derivatives (now properly padded to target_dim)
    derivs <- compute_derivatives(interpolated)
    
    first_deriv_list[[i]] <- derivs$first
    second_deriv_list[[i]] <- derivs$second
  }
  
  # Create time coordinates
  t <- seq(0, 1, length.out = target_dim)
  
  # Convert to WICMAD format (multivariate with 3 channels: original, first deriv, second deriv)
  Y <- lapply(seq_len(n_series), function(i) {
    # Interpolate original series
    original_interp <- interpolate_series(series_data[i, ], target_dim)
    
    # Combine original, first derivative, and second derivative
    matrix(c(original_interp, first_deriv_list[[i]], second_deriv_list[[i]]), 
           nrow = target_dim, ncol = 3)
  })
  
  return(list(
    Y = Y, 
    t = t, 
    labels = labels,
    original_data = series_data,
    first_deriv_data = do.call(rbind, first_deriv_list),
    second_deriv_data = do.call(rbind, second_deriv_list)
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
    warmup = 500,
    pinned_normal = TRUE
  )
  
  # Extract cluster assignments
  cluster_assignments <- wicmad_result$Z[nrow(wicmad_result$Z), ]
  
  # Map WICMAD clusters to binary classification (normal vs anomaly)
  normal_samples <- which(wicmad_data$labels == 0)
  cluster_counts <- table(cluster_assignments[normal_samples])
  normal_cluster <- as.numeric(names(cluster_counts)[which.max(cluster_counts)])
  binary_cluster_assignments <- ifelse(cluster_assignments == normal_cluster, 0, 1)
  
  cat("WICMAD completed. Number of clusters found:", length(unique(cluster_assignments)), "\n")
  cat("Normal cluster ID:", normal_cluster, "\n")
  cat("Binary cluster assignments:\n")
  print(table(binary_cluster_assignments))
  
  # Calculate confusion matrix
  conf_matrix <- table(wicmad_data$labels, binary_cluster_assignments)
  cat("Confusion Matrix:\n")
  print(conf_matrix)
  
  # Calculate accuracy
  accuracy <- sum(diag(conf_matrix)) / sum(conf_matrix)
  cat("Accuracy:", round(accuracy, 4), "\n")
  
  # Calculate ARI
  if (requireNamespace("mclust", quietly = TRUE)) {
    ari <- mclust::adjustedRandIndex(wicmad_data$labels, binary_cluster_assignments)
    cat("Adjusted Rand Index:", round(ari, 4), "\n")
  }
  
  return(list(
    cluster_assignments = cluster_assignments,
    binary_cluster_assignments = binary_cluster_assignments,
    accuracy = accuracy,
    conf_matrix = conf_matrix
  ))
}

# Function to plot overlapped time series
plot_overlapped_data <- function(series_data, labels, title = "Earthquakes - Before Clustering", max_series = NULL) {
  if (is.null(max_series)) {
    max_series <- min(50, nrow(series_data))
  }
  
  # Sample series for plotting
  n_plot <- min(max_series, nrow(series_data))
  plot_indices <- sample(nrow(series_data), n_plot)
  
  # Prepare data for plotting
  plot_data <- data.frame()
  for (i in seq_len(n_plot)) {
    idx <- plot_indices[i]
    series_df <- data.frame(
      time = seq_len(ncol(series_data)),
      value = as.numeric(series_data[idx, ]),
      series_id = paste("Series", idx),
      label = ifelse(labels[idx] == 0, "Normal", "Earthquake")
    )
    plot_data <- rbind(plot_data, series_df)
  }
  
  # Create plot
  p <- ggplot(plot_data, aes(x = .data$time, y = .data$value, color = .data$label, group = .data$series_id)) +
    geom_line(alpha = 0.7, size = 0.5) +
    labs(title = title, x = "Time", y = "Value", color = "Class") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(p)
}

# Function to plot clustering results
plot_clustering_results <- function(series_data, true_labels, cluster_assignments, title = "Earthquakes - After Clustering", max_series = NULL) {
  if (is.null(max_series)) {
    max_series <- min(50, nrow(series_data))
  }
  
  # Sample series for plotting
  n_plot <- min(max_series, nrow(series_data))
  plot_indices <- sample(nrow(series_data), n_plot)
  
  # Prepare data for plotting
  plot_data <- data.frame()
  for (i in seq_len(n_plot)) {
    idx <- plot_indices[i]
    series_df <- data.frame(
      time = seq_len(ncol(series_data)),
      value = as.numeric(series_data[idx, ]),
      series_id = paste("Series", idx),
      cluster = cluster_assignments[idx],
      cluster_name = ifelse(cluster_assignments[idx] == 0, "Normal", "Anomaly")
    )
    plot_data <- rbind(plot_data, series_df)
  }
  
  # Create plot
  p <- ggplot(plot_data, aes(x = .data$time, y = .data$value, color = .data$cluster_name, group = .data$series_id)) +
    geom_line(alpha = 0.7, size = 0.5) +
    labs(title = title, x = "Time", y = "Value", color = "Class") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(p)
}

# Function to plot derivative clustering results
plot_derivative_clustering_results <- function(original_data, first_deriv_data, second_deriv_data, true_labels, cluster_assignments, title = "Earthquakes - After Clustering") {
  # Sample for plotting
  n_plot <- min(30, nrow(original_data))
  plot_indices <- sample(nrow(original_data), n_plot)
  
  # Prepare data for plotting
  plot_data <- data.frame()
  for (i in seq_len(n_plot)) {
    idx <- plot_indices[i]
    
    # Original series
    orig_df <- data.frame(
      time = seq_len(ncol(original_data)),
      value = as.numeric(original_data[idx, ]),
      series_id = paste("Series", idx),
      type = "Original",
      cluster = cluster_assignments[idx],
      cluster_name = ifelse(cluster_assignments[idx] == 0, "Normal", "Anomaly")
    )
    
    # First derivative
    first_df <- data.frame(
      time = seq_len(ncol(first_deriv_data)),
      value = as.numeric(first_deriv_data[idx, ]),
      series_id = paste("Series", idx),
      type = "First Derivative",
      cluster = cluster_assignments[idx],
      cluster_name = ifelse(cluster_assignments[idx] == 0, "Normal", "Anomaly")
    )
    
    # Second derivative
    second_df <- data.frame(
      time = seq_len(ncol(second_deriv_data)),
      value = as.numeric(second_deriv_data[idx, ]),
      series_id = paste("Series", idx),
      type = "Second Derivative",
      cluster = cluster_assignments[idx],
      cluster_name = ifelse(cluster_assignments[idx] == 0, "Normal", "Anomaly")
    )
    
    plot_data <- rbind(plot_data, orig_df, first_df, second_df)
  }
  
  # Create plot
  p <- ggplot(plot_data, aes(x = .data$time, y = .data$value, color = .data$cluster_name, group = .data$series_id)) +
    geom_line(alpha = 0.7, size = 0.5) +
    facet_wrap(~ type, scales = "free_y", ncol = 1) +
    labs(title = title, x = "Time", y = "Value", color = "Class") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(p)
}


# Main analysis function
main <- function() {
  cat("=== Earthquakes Dataset Analysis ===\n")
  
  # Set data directory
  data_dir <- "../data/Earthquakes"
  
  # Load data
  cat("\n1. Loading earthquakes dataset...\n")
  data <- load_earthquakes_data(data_dir, imbalance_ratio = 0.05)
  
  # Limit to 30 observations for testing
  if (nrow(data$train_data) > 30) {
    cat("Limiting to 30 observations for testing...\n")
    indices <- sample(nrow(data$train_data), 30)
    data$train_data <- data$train_data[indices, ]
    data$train_labels <- data$train_labels[indices]
  }
  
  cat("Label distribution:\n")
  print(table(data$train_labels))
  cat("Majority class used as normal:", data$majority_class, "\n")
  cat("First few labels:", head(data$train_labels, 10), "\n")
  
  # Plot original data (overlapped)
  cat("\n2. Creating original data visualization...\n")
  original_plot <- plot_overlapped_data(data$train_data, data$train_labels, 
                                       "Earthquakes Dataset - Before Clustering")
  
  # Save original data plot
  ggsave("../plots/earthquakes/earthquakes_imbalanced_original_data.pdf", original_plot, 
         width = 10, height = 6, dpi = 300)
  cat("Saved: ../plots/earthquakes/earthquakes_imbalanced_original_data.pdf\n")
  
  # Prepare data for different analyses
  cat("\n3. Preparing data for WICMAD analysis...\n")
  
  # Raw time series analysis
  raw_wicmad_data <- prepare_wicmad_data(data$train_data, data$train_labels, data$target_dim)
  
  # Derivative analysis
  deriv_data_prep <- prepare_derivative_data(data$train_data, data$train_labels, data$target_dim)
  
  
  # Run WICMAD analyses
  cat("\n4. Running WICMAD clustering...\n")
  
  # Raw time series clustering
  raw_results <- run_wicmad_analysis(raw_wicmad_data, "Raw Time Series")
  
  # Derivative clustering
  deriv_results <- run_wicmad_analysis(deriv_data_prep, "Derivatives")
  
  
  # Create clustering visualizations
  cat("\n5. Creating clustering visualizations...\n")
  
  # Raw time series clustering plot
  raw_clustering_plot <- plot_clustering_results(
    data$train_data,
    data$train_labels,
    raw_results$cluster_assignments,
    "Earthquakes Dataset - After Clustering"
  )
  
  # Derivative clustering plot
  deriv_clustering_plot <- plot_derivative_clustering_results(
    deriv_data_prep$original_data,
    deriv_data_prep$first_deriv_data,
    deriv_data_prep$second_deriv_data,
    data$train_labels,
    deriv_results$cluster_assignments,
    "Earthquakes Dataset - After Clustering"
  )
  
  
  # Save clustering plots
  ggsave("../plots/earthquakes/earthquakes_raw_curves_clustering.pdf", raw_clustering_plot, 
         width = 10, height = 6, dpi = 300)
  cat("Saved: ../plots/earthquakes/earthquakes_raw_curves_clustering.pdf\n")
  
  ggsave("../plots/earthquakes/earthquakes_derivatives_clustering.pdf", deriv_clustering_plot, 
         width = 10, height = 8, dpi = 300)
  cat("Saved: ../plots/earthquakes/earthquakes_derivatives_clustering.pdf\n")
  
  
  cat("\n=== Analysis Complete ===\n")
  cat("Results summary:\n")
  cat("- Raw time series accuracy:", round(raw_results$accuracy, 4), "\n")
  cat("- Derivatives accuracy:", round(deriv_results$accuracy, 4), "\n")
}

# Run the analysis
main()
