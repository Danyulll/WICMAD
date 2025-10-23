#!/usr/bin/env Rscript

# Tecator Dataset Analysis with WICMAD
# This script loads the Tecator dataset, runs WICMAD clustering, and visualizes results

# Set working directory to scripts folder
# setwd("scripts")

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(cluster)
library(fda)  # For FPCA
library(mclust)  # For adjustedRandIndex
library(numDeriv)

# Load WICMAD package
library(devtools)
devtools::load_all()

# Function to load Tecator dataset and create imbalanced data
load_tecator_data <- function(imbalance_ratio = 0.05) {
  cat("Loading Tecator dataset...\n")
  
  # Load the Tecator data
  data(Tecator, package = "fsemipar")
  
  # Extract the data components
  fat_content <- Tecator$fat
  protein_content <- Tecator$protein
  moisture_content <- Tecator$moisture
  absor_spectra <- Tecator$absor.spectra
  absor_spectra1 <- Tecator$absor.spectra1
  absor_spectra2 <- Tecator$absor.spectra2
  
  cat("Dataset dimensions:\n")
  cat("Number of samples:", length(fat_content), "\n")
  cat("Absorbance spectra dimensions:", dim(absor_spectra), "\n")
  cat("Fat content range:", range(fat_content), "\n")
  cat("Protein content range:", range(protein_content), "\n")
  cat("Moisture content range:", range(moisture_content), "\n")
  
  # Set target dimension to 16 for testing
  target_dim <- 16
  cat("Using target dimension for testing:", target_dim, "\n")
  
  # Create binary labels based on fat content (high fat vs low fat)
  # Use median as threshold for binary classification
  fat_threshold <- median(fat_content)
  labels <- ifelse(fat_content > fat_threshold, 1, 0)  # 1 = high fat, 0 = low fat
  
  cat("Original distribution:\n")
  cat("Low fat samples:", sum(labels == 0), "\n")
  cat("High fat samples:", sum(labels == 1), "\n")
  cat("Fat threshold:", round(fat_threshold, 2), "\n")
  
  # Create imbalanced dataset with 5% anomalies
  cat("Creating imbalanced dataset with 5% anomalies...\n")
  
  # Identify majority and minority classes
  class_counts <- table(labels)
  majority_class <- as.numeric(names(class_counts)[which.max(class_counts)])
  minority_classes <- as.numeric(names(class_counts)[names(class_counts) != majority_class])
  
  cat("Majority class:", majority_class, "with", max(class_counts), "samples\n")
  cat("Minority classes:", paste(minority_classes, collapse = ", "), "\n")
  
  # Create binary labels: majority = 0 (normal), minority = 1 (anomaly)
  binary_labels <- ifelse(labels == majority_class, 0, 1)
  
  # Calculate how many anomalies we need for 5%
  total_samples <- length(binary_labels)
  n_anomalies_needed <- round(total_samples * imbalance_ratio)
  n_normal_needed <- total_samples - n_anomalies_needed
  
  cat("Target: Normal samples:", n_normal_needed, "Anomaly samples:", n_anomalies_needed, "\n")
  
  # Get indices for each class
  normal_indices <- which(binary_labels == 0)
  anomaly_indices <- which(binary_labels == 1)
  
  # Sample the required number of each class
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
  imbalanced_labels <- binary_labels[selected_indices]
  imbalanced_spectra <- absor_spectra[selected_indices, ]
  imbalanced_spectra1 <- absor_spectra1[selected_indices, ]
  imbalanced_spectra2 <- absor_spectra2[selected_indices, ]
  imbalanced_fat <- fat_content[selected_indices]
  imbalanced_protein <- protein_content[selected_indices]
  imbalanced_moisture <- moisture_content[selected_indices]
  
  cat("Final imbalanced dataset:\n")
  cat("Total samples:", length(imbalanced_labels), "\n")
  cat("Normal (majority class):", sum(imbalanced_labels == 0), "\n")
  cat("Anomaly (minority class):", sum(imbalanced_labels == 1), "\n")
  cat("Anomaly percentage:", round(mean(imbalanced_labels == 1) * 100, 1), "%\n")
  
  return(list(
    train_labels = imbalanced_labels,
    train_spectra = imbalanced_spectra,
    train_spectra1 = imbalanced_spectra1,
    train_spectra2 = imbalanced_spectra2,
    train_fat = imbalanced_fat,
    train_protein = imbalanced_protein,
    train_moisture = imbalanced_moisture,
    test_labels = NULL,
    test_spectra = NULL,
    target_dim = target_dim
  ))
}

# Function to interpolate time series to target dimension
interpolate_series <- function(series, target_dim = 16) {
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

# Function to perform functional PCA and select components for 90% variance
perform_fpca <- function(spectra, target_dim = 16) {
  cat("Performing functional PCA analysis...\n")
  
  # Interpolate all spectra to same length
  n_series <- nrow(spectra)
  interpolated_spectra <- matrix(0, nrow = n_series, ncol = target_dim)
  
  for (i in seq_len(n_series)) {
    interpolated_spectra[i, ] <- interpolate_series(spectra[i, ], target_dim)
  }
  
  # Create time points
  t <- seq(0, 1, length.out = target_dim)
  
  # Perform functional PCA using fda package
  # Create functional data object
  fd_obj <- fda::Data2fd(argvals = t, y = t(interpolated_spectra))
  
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

# Function to prepare data for WICMAD (raw spectra)
prepare_wicmad_data <- function(spectra, labels, target_dim = 64) {
  cat("Preparing raw spectra data for WICMAD...\n")
  
  # Interpolate each spectrum to target dimensions
  n_series <- nrow(spectra)
  interpolated_spectra <- matrix(0, nrow = n_series, ncol = target_dim)
  
  for (i in seq_len(n_series)) {
    interpolated_spectra[i, ] <- interpolate_series(spectra[i, ], target_dim)
  }
  
  # Create time coordinates (normalized to [0,1])
  t <- seq(0, 1, length.out = target_dim)
  
  # Convert to list format for WICMAD (each element is a PxM matrix)
  # For univariate data, we need P=target_dim (time points) and M=1 (channels)
  Y <- lapply(seq_len(n_series), function(i) {
    matrix(interpolated_spectra[i, ], nrow = target_dim, ncol = 1)
  })
  
  return(list(Y = Y, t = t, labels = labels))
}

# Function to prepare derivative data for WICMAD
prepare_derivative_data <- function(original_spectra, spectra1, spectra2, labels) {
  cat("Preparing derivative data for WICMAD...\n")
  
  n_series <- nrow(original_spectra)
  
  # Interpolate original data to 16 dimensions first
  original_interp <- matrix(0, nrow = n_series, ncol = 16)
  for (i in seq_len(n_series)) {
    original_interp[i, ] <- interpolate_series(original_spectra[i, ], 16)
  }
  
  # Interpolate first derivative
  first_deriv_interp <- matrix(0, nrow = n_series, ncol = 16)
  for (i in seq_len(n_series)) {
    first_deriv_interp[i, ] <- interpolate_series(spectra1[i, ], 16)
  }
  
  # Interpolate second derivative
  second_deriv_interp <- matrix(0, nrow = n_series, ncol = 16)
  for (i in seq_len(n_series)) {
    second_deriv_interp[i, ] <- interpolate_series(spectra2[i, ], 16)
  }
  
  # Create time coordinates
  t <- seq(0, 1, length.out = 16)
  
  # Convert to list format for WICMAD (each element is a PxM matrix)
  # P=16 (time points), M=3 (original + first + second derivatives)
  Y <- lapply(seq_len(n_series), function(i) {
    matrix(c(original_interp[i, ], first_deriv_interp[i, ], second_deriv_interp[i, ]), nrow = 16, ncol = 3)
  })
  
  return(list(Y = Y, t = t, labels = labels, 
              original_data = original_interp,
              first_deriv_data = first_deriv_interp,
              second_deriv_data = second_deriv_interp))
}

# Function to prepare FPCA data for WICMAD
prepare_fpca_data <- function(spectra, labels) {
  cat("Preparing FPCA data for WICMAD...\n")
  
  # Perform functional PCA
  fpca_result <- perform_fpca(spectra, 16)
  
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

# Function to plot clustering results
plot_clustering_results <- function(spectra, true_labels, cluster_assignments, title = "Tecator Dataset", max_series = NULL) {
  cat("Creating clustering results plot...\n")
  
  # Use all series if max_series is NULL
  if (is.null(max_series)) {
    series_indices <- seq_len(nrow(spectra))
  } else {
    n_series <- min(max_series, nrow(spectra))
    series_indices <- sample(seq_len(nrow(spectra)), n_series)
  }
  
  # Create plotting data
  plot_data <- data.frame()
  
  for (i in series_indices) {
    series_data <- spectra[i, ]
    time_points <- seq(0, 1, length.out = length(series_data))
    
    series_df <- data.frame(
      time = time_points,
      value = series_data,
      series_id = i,
      true_label = true_labels[i],
      cluster = cluster_assignments[i],
      cluster_name = paste("Cluster", cluster_assignments[i])
    )
    
    plot_data <- rbind(plot_data, series_df)
  }
  
  # Create overlapped plot
  p <- ggplot(plot_data, aes(x = .data$time, y = .data$value, color = .data$cluster_name, group = .data$series_id)) +
    geom_line(alpha = 0.6, linewidth = 0.5) +
    labs(title = title,
         x = "Wavelength (normalized)", y = "Absorbance", color = "Cluster") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(p)
}

# Function to plot overlapped data
plot_overlapped_data <- function(spectra, labels, title = "Tecator Dataset", max_series = NULL) {
  cat("Creating overlapped data plot...\n")
  
  # Use all series if max_series is NULL
  if (is.null(max_series)) {
    series_indices <- seq_len(nrow(spectra))
  } else {
    n_series <- min(max_series, nrow(spectra))
    series_indices <- sample(seq_len(nrow(spectra)), n_series)
  }
  
  # Create plotting data
  plot_data <- data.frame()
  
  for (i in series_indices) {
    series_data <- spectra[i, ]
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
         x = "Wavelength (normalized)", y = "Absorbance", color = "Class") +
    scale_color_manual(values = c("Normal" = "blue", "Anomaly" = "red")) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(p)
}

# Function to plot derivative clustering results (3 separate channels)
plot_derivative_clustering_results <- function(original_data, first_deriv_data, second_deriv_data, true_labels, cluster_assignments, title = "Tecator - After Clustering", max_series = NULL) {
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
      cluster_name = ifelse(cluster_assignments[i] == 0, "Normal", "Anomaly")
    )
    
    # First derivative
    first_df <- data.frame(
      time = time_points,
      value = first_deriv_data[i, ],
      channel = "1st Derivative",
      series_id = i,
      true_label = true_labels[i],
      cluster = cluster_assignments[i],
      cluster_name = ifelse(cluster_assignments[i] == 0, "Normal", "Anomaly")
    )
    
    # Second derivative
    second_df <- data.frame(
      time = time_points,
      value = second_deriv_data[i, ],
      channel = "2nd Derivative",
      series_id = i,
      true_label = true_labels[i],
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
         x = "Wavelength (normalized)", y = "Value", color = "Class") +
    scale_color_manual(values = c("Normal" = "blue", "Anomaly" = "red")) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(p)
}

# Function to plot FPCA clustering results (showing functional components)
plot_fpca_clustering_results <- function(fpca_result, true_labels, cluster_assignments, title = "Tecator - After Clustering", max_series = NULL) {
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
        cluster_name = ifelse(cluster_assignments[i] == 0, "Normal", "Anomaly"),
        stringsAsFactors = FALSE
      )
      
      plot_data <- rbind(plot_data, component_df)
    }
  }
  
  # Create faceted plot for functional components
  p <- ggplot(plot_data, aes(x = .data$time, y = .data$value, color = .data$cluster_name, group = .data$series_id)) +
    geom_line(alpha = 0.6, linewidth = 0.5) +
    facet_wrap(~ .data$component_name, scales = "free_y", ncol = 2) +
    labs(title = title,
         x = "Wavelength (normalized)", y = "Component Value", color = "Class") +
    scale_color_manual(values = c("Normal" = "blue", "Anomaly" = "red")) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(p)
}

# Function to run WICMAD analysis on a dataset
run_wicmad_analysis <- function(data_prep, data_name, spectra, labels, spectra1 = NULL, spectra2 = NULL) {
  cat(paste("\n=== Running WICMAD on", data_name, "===\n"))
  
  # Prepare data
  if (data_name == "Derivatives") {
    wicmad_data <- data_prep(spectra, spectra1, spectra2, labels)
  } else {
    wicmad_data <- data_prep(spectra, labels)
  }
  
  cat("Number of series:", length(wicmad_data$Y), "\n")
  cat("Time points:", length(wicmad_data$t), "\n")
  cat("Series dimensions:", nrow(wicmad_data$Y[[1]]), "x", ncol(wicmad_data$Y[[1]]), "\n")
  
  # Run WICMAD
  cat("Running WICMAD clustering...\n")
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
  metrics <- calculate_clustering_metrics(wicmad_data$labels, binary_cluster_assignments)
  
  cat("Clustering Performance:\n")
  cat("Macro Precision:", round(metrics$Macro_Precision, 4), "\n")
  cat("Macro Recall:", round(metrics$Macro_Recall, 4), "\n")
  cat("Macro Specificity:", round(metrics$Macro_Specificity, 4), "\n")
  cat("Macro F1 Score:", round(metrics$Macro_F1, 4), "\n")
  cat("Overall Accuracy:", round(metrics$Overall_Accuracy, 4), "\n")
  cat("Adjusted Rand Index:", round(metrics$ARI, 4), "\n")
  
  return(list(
    metrics = metrics,
    cluster_assignments = binary_cluster_assignments,
    wicmad_result = wicmad_result,
    wicmad_data = wicmad_data
  ))
}

# Main analysis
main <- function() {
  cat("=== Tecator Dataset Analysis with WICMAD ===\n\n")
  
  # Load data with 5% high fat class
  cat("1. Loading Tecator dataset with 5% high fat class...\n")
  data <- load_tecator_data(imbalance_ratio = 0.05)
  
  # Get the calculated target dimension
  target_dim <- data$target_dim
  cat("Using target dimension P =", target_dim, "(highest power of 2 available)\n")
  
  # Limit to 30 observations for testing
  if (nrow(data$train_spectra) > 30) {
    cat("Limiting to 30 observations for testing...\n")
    indices <- sample(nrow(data$train_spectra), 30)
    data$train_spectra <- data$train_spectra[indices, ]
    data$train_spectra1 <- data$train_spectra1[indices, ]
    data$train_spectra2 <- data$train_spectra2[indices, ]
    data$train_labels <- data$train_labels[indices]
    data$train_fat <- data$train_fat[indices]
    data$train_protein <- data$train_protein[indices]
    data$train_moisture <- data$train_moisture[indices]
  }
  
  # Print data summary
  cat("Training data shape:", nrow(data$train_spectra), "x", ncol(data$train_spectra), "\n")
  cat("Label distribution:\n")
  print(table(data$train_labels))
  cat("First few labels:", head(data$train_labels, 10), "\n")
  
  # Plot original data (overlapped)
  cat("\n2. Creating original data visualization...\n")
  original_plot <- plot_overlapped_data(data$train_spectra, data$train_labels, 
                                       "Tecator Dataset - Before Clustering")
  
  # Save original data plot
  pdf("../plots/tecator/tecator_imbalanced_original_data.pdf", width = 12, height = 8)
  print(original_plot)
  dev.off()
  cat("Imbalanced original data plot saved to ../plots/tecator/tecator_imbalanced_original_data.pdf\n")
  
  # Run analysis on raw spectra
  cat("\n3. Analyzing raw spectra data...\n")
  raw_results <- run_wicmad_analysis(
    function(spectra, labels) prepare_wicmad_data(spectra, labels, target_dim), 
    "Raw Spectra", 
    data$train_spectra, 
    data$train_labels
  )
  
  # Plot raw spectra clustering results
  raw_clustering_plot <- plot_clustering_results(
    data$train_spectra, 
    data$train_labels, 
    raw_results$cluster_assignments,
    "Tecator Dataset - After Clustering"
  )
  
  pdf("../plots/tecator/tecator_raw_spectra_clustering.pdf", width = 12, height = 8)
  print(raw_clustering_plot)
  dev.off()
  cat("Raw spectra clustering plot saved to ../plots/tecator/tecator_raw_spectra_clustering.pdf\n")
  
  # Run analysis on derivatives
  cat("\n4. Analyzing derivative data...\n")
  deriv_results <- run_wicmad_analysis(
    prepare_derivative_data, 
    "Derivatives", 
    data$train_spectra, 
    data$train_labels,
    data$train_spectra1,
    data$train_spectra2
  )
  
  # Plot derivative clustering results (3 separate channels)
  deriv_data_prep <- prepare_derivative_data(data$train_spectra, data$train_spectra1, data$train_spectra2, data$train_labels)
  deriv_clustering_plot <- plot_derivative_clustering_results(
    deriv_data_prep$original_data,
    deriv_data_prep$first_deriv_data,
    deriv_data_prep$second_deriv_data,
    deriv_results$cluster_assignments, 
    deriv_results$cluster_assignments,
    "Tecator Dataset - After Clustering"
  )
  
  pdf("../plots/tecator/tecator_derivatives_clustering.pdf", width = 12, height = 8)
  print(deriv_clustering_plot)
  dev.off()
  cat("Derivatives clustering plot saved to ../plots/tecator/tecator_derivatives_clustering.pdf\n")
  
  # Run analysis on FPCA data
  cat("\n5. Analyzing FPCA data...\n")
  fpca_results <- run_wicmad_analysis(
    prepare_fpca_data, 
    "FPCA", 
    data$train_spectra, 
    data$train_labels
  )
  
  # Plot FPCA clustering results (showing functional components)
  # Get FPCA data from the preparation step
  fpca_data_prep <- prepare_fpca_data(data$train_spectra, data$train_labels)
  fpca_clustering_plot <- plot_fpca_clustering_results(
    fpca_data_prep$fpca_result, 
    fpca_results$cluster_assignments, 
    fpca_results$cluster_assignments,
    "Tecator Dataset - After Clustering"
  )
  
  pdf("../plots/tecator/tecator_fpca_clustering.pdf", width = 12, height = 8)
  print(fpca_clustering_plot)
  dev.off()
  cat("FPCA clustering plot saved to ../plots/tecator/tecator_fpca_clustering.pdf\n")
  
  # Print final summary
  cat("\n=== Imbalanced Analysis Complete ===\n")
  cat("Generated files:\n")
  cat("- ../plots/tecator/tecator_imbalanced_original_data.pdf: Imbalanced spectra (5% high fat, overlapped)\n")
  cat("- ../plots/tecator/tecator_raw_spectra_clustering.pdf: Raw spectra clustering results\n")
  cat("- ../plots/tecator/tecator_derivatives_clustering.pdf: Derivatives clustering results\n")
  cat("- ../plots/tecator/tecator_fpca_clustering.pdf: FPCA clustering results\n")
  
  cat("\nPerformance Comparison:\n")
  cat("Raw Spectra - Precision:", round(raw_results$metrics$Macro_Precision, 4), 
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
