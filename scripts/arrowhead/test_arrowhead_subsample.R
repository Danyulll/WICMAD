#!/usr/bin/env Rscript

# Test script for ArrowHead subsample with wavelet C++ implementation
# This script runs WICMAD on a small subsample to verify wavelet code works correctly

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(cluster)
library(mclust)
library(numDeriv)

# Load WICMAD package
library(devtools)
devtools::load_all()

# Function to load ArrowHead dataset (simplified for testing)
load_arrowhead_subsample <- function(data_dir, n_samples = 20) {
  cat("Loading ArrowHead dataset subsample...\n")
  
  # Load training data
  train_file <- file.path(data_dir, "ArrowHead_TRAIN.ts")
  
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
  cat("Total series available:", n_series, "\n")
  
  # Take a subsample
  if (n_samples > n_series) {
    n_samples <- n_series
    cat("Requested sample size larger than available data, using all", n_samples, "series\n")
  }
  
  # Sample indices
  sample_indices <- sample(seq_len(n_series), n_samples)
  cat("Using subsample of", n_samples, "series\n")
  
  # Parse selected series
  series_list <- list()
  labels <- integer(n_samples)
  
  for (i in seq_along(sample_indices)) {
    line_idx <- sample_indices[i]
    line <- data_lines[line_idx]
    
    # Split by colon to separate data from label
    parts <- strsplit(line, ":")[[1]]
    if (length(parts) != 2) {
      stop("Invalid data format in line ", line_idx, " - expected 2 parts, got ", length(parts))
    }
    
    # Parse the data part (comma-separated values)
    data_values <- as.numeric(strsplit(parts[1], ",")[[1]])
    
    # Parse the label part (integer)
    labels[i] <- as.integer(parts[2])
    
    # Store as a single column matrix (univariate time series)
    series_list[[i]] <- matrix(data_values, ncol = 1)
  }
  
  # Find the majority class to use as normal
  class_counts <- table(labels)
  majority_class <- as.numeric(names(class_counts)[which.max(class_counts)])
  minority_classes <- as.numeric(names(class_counts)[class_counts != max(class_counts)])
  
  cat("Subsample class distribution:\n")
  print(class_counts)
  cat("Majority class (normal):", majority_class, "\n")
  cat("Minority classes (anomalies):", paste(minority_classes, collapse = ", "), "\n")
  
  # Convert to anomaly detection format (majority class = normal, others = anomaly)
  anomaly_labels <- ifelse(labels == majority_class, 0, 1)
  
  cat("Final subsample dataset:\n")
  cat("Total samples:", length(series_list), "\n")
  cat("Normal (Class", majority_class, "):", sum(anomaly_labels == 0), "\n")
  cat("Anomaly (Classes", paste(minority_classes, collapse = ", "), "):", sum(anomaly_labels == 1), "\n")
  cat("Anomaly percentage:", round(mean(anomaly_labels == 1) * 100, 1), "%\n")
  
  return(list(
    train_series = series_list,
    train_labels = anomaly_labels
  ))
}

# Function to interpolate time series to target dimension
interpolate_series <- function(series, target_dim = 32) {
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
  
  # Use p=32 for higher resolution
  p_dim <- 32
  cat("Using p parameter:", p_dim, "\n")
  
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

# Function to test wavelet functions
test_wavelet_functions <- function(series_list) {
  cat("\n=== Testing Wavelet Functions ===\n")
  
  # Test on first series
  test_series <- series_list[[1]][, 1]  # Extract univariate data
  cat("Testing wavelet functions on series of length:", length(test_series), "\n")
  
  # Test forward transform
  cat("\n1. Testing forward wavelet transform...\n")
  start_time <- Sys.time()
  tryCatch({
    wt_result <- wt_forward_1d(test_series, wf = "la8", J = 4, boundary = "periodic")
    forward_time <- as.numeric(Sys.time() - start_time, units = "secs")
    cat("✓ Forward transform successful\n")
    cat("  Time:", round(forward_time, 4), "seconds\n")
    cat("  Coefficients length:", length(wt_result$coeff), "\n")
    cat("  J levels:", wt_result$map$J, "\n")
  }, error = function(e) {
    cat("✗ Forward transform failed:", e$message, "\n")
    return(FALSE)
  })
  
  # Test inverse transform
  cat("\n2. Testing inverse wavelet transform...\n")
  start_time <- Sys.time()
  tryCatch({
    reconstructed <- wt_inverse_1d(wt_result$coeff, wt_result$map)
    inverse_time <- as.numeric(Sys.time() - start_time, units = "secs")
    cat("✓ Inverse transform successful\n")
    cat("  Time:", round(inverse_time, 4), "seconds\n")
    cat("  Reconstructed length:", length(reconstructed), "\n")
    
    # Check reconstruction accuracy
    mse <- mean((test_series - reconstructed)^2)
    cat("  Reconstruction MSE:", round(mse, 8), "\n")
    if (mse < 1e-10) {
      cat("  ✓ Perfect reconstruction!\n")
    } else if (mse < 1e-6) {
      cat("  ✓ Excellent reconstruction (numerical precision)\n")
    } else {
      cat("  ⚠ Some reconstruction error\n")
    }
  }, error = function(e) {
    cat("✗ Inverse transform failed:", e$message, "\n")
    return(FALSE)
  })
  
  # Test C++ implementation availability
  cat("\n3. Checking C++ implementation...\n")
  if (exists("fast_wavelet_forward_1d")) {
    cat("✓ C++ function 'fast_wavelet_forward_1d' is available\n")
    
    # Test C++ implementation
    start_time <- Sys.time()
    tryCatch({
      cpp_result <- fast_wavelet_forward_1d(test_series, wf = "la8", J = 4, boundary = "periodic")
      cpp_time <- as.numeric(Sys.time() - start_time, units = "secs")
      cat("✓ C++ forward transform successful\n")
      cat("  Time:", round(cpp_time, 4), "seconds\n")
      
      # Compare with R implementation
      coeff_diff <- mean(abs(wt_result$coeff - cpp_result$coeff))
      cat("  Coefficient difference from R:", round(coeff_diff, 8), "\n")
      
      if (coeff_diff < 1e-6) {
        cat("  ✓ C++ and R results match!\n")
      } else {
        cat("  ⚠ Some difference between C++ and R results\n")
      }
      
      # Performance comparison
      if (forward_time > 0 && cpp_time > 0) {
        speedup <- forward_time / cpp_time
        cat("  C++ speedup:", round(speedup, 2), "x\n")
      }
      
    }, error = function(e) {
      cat("✗ C++ implementation failed:", e$message, "\n")
      cat("  Falling back to R implementation\n")
    })
  } else {
    cat("✗ C++ function 'fast_wavelet_forward_1d' is NOT available\n")
    cat("  Using R implementation only\n")
  }
  
  return(TRUE)
}

# Main test function
main <- function() {
  cat("=== ArrowHead Subsample Test with Wavelet C++ Implementation ===\n\n")
  
  # Set data directory
  data_dir <- "../../data/ArrowHead"
  
  # Load subsample data
  cat("1. Loading ArrowHead subsample...\n")
  data <- load_arrowhead_subsample(data_dir, n_samples = 20)
  
  # Test wavelet functions
  test_wavelet_functions(data$train_series)
  
  # Prepare data for WICMAD
  cat("\n2. Preparing data for WICMAD...\n")
  wicmad_data <- prepare_wicmad_data(data$train_series, data$train_labels)
  
  cat("Number of series:", length(wicmad_data$Y), "\n")
  cat("Time points:", length(wicmad_data$t), "\n")
  cat("Series dimensions:", nrow(wicmad_data$Y[[1]]), "x", ncol(wicmad_data$Y[[1]]), "\n")
  
  # Run WICMAD clustering with minimal iterations for testing
  cat("\n3. Running WICMAD clustering (test mode - 50 iterations)...\n")
  start_time <- Sys.time()
  
  tryCatch({
    wicmad_result <- wicmad(
      Y = wicmad_data$Y,
      t = wicmad_data$t,
      n_iter = 50,   # Minimal for testing
      burn = 25,     # Minimal for testing
      thin = 1,
      warmup_iters = 5,  # Minimal for testing
      unpin = FALSE,
      alpha_prior = c(5, 1),  # Adjusted from c(50, 1)
      kappa_pi = 0.3,
      c2 = 0.5,
      tau_pi = 10,
      diagnostics = TRUE
    )
    
    total_time <- as.numeric(Sys.time() - start_time, units = "secs")
    
    cat("✓ WICMAD completed successfully\n")
    cat("  Total time:", round(total_time, 2), "seconds\n")
    cat("  Number of iterations:", 50, "\n")
    cat("  Time per iteration:", round(total_time / 50, 4), "seconds\n")
    
    # Check number of clusters
    if (!is.null(wicmad_result$K_occ)) {
      cat("  Number of occupied clusters (last 10 iterations):\n")
      last_10 <- tail(wicmad_result$K_occ, 10)
      cat("    Range:", min(last_10), "-", max(last_10), "\n")
      cat("    Mean:", round(mean(last_10), 1), "\n")
      cat("    Median:", median(last_10), "\n")
    }
    
    # Check if C++ implementation was used
    cat("\n4. Checking implementation usage...\n")
    if (exists("fast_icm_loglik_curve_eigen")) {
      cat("✓ C++ function 'fast_icm_loglik_curve_eigen' is available\n")
    } else {
      cat("✗ C++ function 'fast_icm_loglik_curve_eigen' is NOT available\n")
    }
    
    if (exists("fast_wavelet_forward_1d")) {
      cat("✓ C++ function 'fast_wavelet_forward_1d' is available\n")
    } else {
      cat("✗ C++ function 'fast_wavelet_forward_1d' is NOT available\n")
    }
    
    # Performance analysis
    cat("\n=== Performance Analysis ===\n")
    cat("With alpha_prior = c(5, 1):\n")
    cat("- Expected to produce fewer clusters than c(50, 1)\n")
    cat("- Should be faster due to fewer cluster computations\n")
    cat("- C++ implementation should provide speedup\n")
    
    # Test wavelet performance on multiple series
    cat("\n5. Testing wavelet performance on multiple series...\n")
    n_test_series <- min(5, length(data$train_series))
    test_indices <- sample(seq_len(length(data$train_series)), n_test_series)
    
    total_wavelet_time <- 0
    for (i in test_indices) {
      series_data <- data$train_series[[i]][, 1]
      start_time <- Sys.time()
      wt_result <- wt_forward_1d(series_data, wf = "la8", J = 4, boundary = "periodic")
      reconstructed <- wt_inverse_1d(wt_result$coeff, wt_result$map)
      wavelet_time <- as.numeric(Sys.time() - start_time, units = "secs")
      total_wavelet_time <- total_wavelet_time + wavelet_time
    }
    
    cat("  Tested", n_test_series, "series\n")
    cat("  Total wavelet time:", round(total_wavelet_time, 4), "seconds\n")
    cat("  Average per series:", round(total_wavelet_time / n_test_series, 4), "seconds\n")
    
  }, error = function(e) {
    cat("✗ WICMAD failed:", e$message, "\n")
    cat("This might be due to compilation issues or other errors.\n")
  })
  
  cat("\n=== Test Complete ===\n")
  cat("Summary:\n")
  cat("- Wavelet functions tested for correctness\n")
  cat("- C++ implementation availability checked\n")
  cat("- WICMAD algorithm tested on subsample\n")
  cat("- Performance metrics collected\n")
}

# Run the test
main()
