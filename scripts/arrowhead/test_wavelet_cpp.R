#!/usr/bin/env Rscript

# Test script for C++ wavelet implementation
# This script tests the new C++ wavelet functions for correctness and performance

# Load required libraries
library(devtools)
devtools::load_all()

# Test function for C++ wavelet implementation
test_wavelet_cpp <- function() {
  cat("=== Testing C++ Wavelet Implementation ===\n\n")
  
  # Test data
  set.seed(123)
  P <- 32  # Test with p=32 as we're using now
  y <- rnorm(P)
  
  cat("Test data created:\n")
  cat("- Signal length:", P, "\n")
  cat("- Signal range:", round(range(y), 4), "\n\n")
  
  # Test if C++ functions are available
  cat("Testing C++ function availability...\n")
  if (exists("fast_wavelet_forward_1d")) {
    cat("✓ C++ function 'fast_wavelet_forward_1d' is available\n")
  } else {
    cat("✗ C++ function 'fast_wavelet_forward_1d' is NOT available\n")
    cat("  This means the C++ code hasn't been compiled yet.\n")
    cat("  You need to run: Rcpp::compileAttributes() and then rebuild the package.\n")
    return(FALSE)
  }
  
  if (exists("fast_wavelet_inverse_1d")) {
    cat("✓ C++ function 'fast_wavelet_inverse_1d' is available\n")
  } else {
    cat("✗ C++ function 'fast_wavelet_inverse_1d' is NOT available\n")
    return(FALSE)
  }
  
  # Test forward transform
  cat("\nTesting forward wavelet transform...\n")
  start_time <- Sys.time()
  tryCatch({
    result_cpp <- fast_wavelet_forward_1d(y, wf = "la8", J = 4, boundary = "periodic")
    cpp_time <- as.numeric(Sys.time() - start_time, units = "secs")
    cat("✓ C++ forward transform successful\n")
    cat("  Time:", round(cpp_time, 4), "seconds\n")
    cat("  Coefficients length:", length(result_cpp$coeff), "\n")
    cat("  J levels:", result_cpp$map$J, "\n")
  }, error = function(e) {
    cat("✗ C++ forward transform failed:", e$message, "\n")
    return(FALSE)
  })
  
  # Test inverse transform
  cat("\nTesting inverse wavelet transform...\n")
  start_time <- Sys.time()
  tryCatch({
    reconstructed <- fast_wavelet_inverse_1d(result_cpp$coeff, result_cpp$map)
    cpp_inv_time <- as.numeric(Sys.time() - start_time, units = "secs")
    cat("✓ C++ inverse transform successful\n")
    cat("  Time:", round(cpp_inv_time, 4), "seconds\n")
    cat("  Reconstructed length:", length(reconstructed), "\n")
    
    # Check reconstruction accuracy
    mse <- mean((y - reconstructed)^2)
    cat("  Reconstruction MSE:", round(mse, 8), "\n")
    if (mse < 1e-10) {
      cat("  ✓ Perfect reconstruction!\n")
    } else {
      cat("  ⚠ Some reconstruction error (may be due to numerical precision)\n")
    }
  }, error = function(e) {
    cat("✗ C++ inverse transform failed:", e$message, "\n")
    return(FALSE)
  })
  
  # Compare with R implementation
  cat("\nComparing with R implementation...\n")
  start_time <- Sys.time()
  tryCatch({
    result_r <- wt_forward_1d(y, wf = "la8", J = 4, boundary = "periodic")
    r_time <- as.numeric(Sys.time() - start_time, units = "secs")
    cat("✓ R implementation successful\n")
    cat("  Time:", round(r_time, 4), "seconds\n")
    
    # Compare results
    coeff_diff <- mean(abs(result_cpp$coeff - result_r$coeff))
    cat("  Coefficient difference:", round(coeff_diff, 8), "\n")
    
    if (coeff_diff < 1e-6) {
      cat("  ✓ C++ and R results match!\n")
    } else {
      cat("  ⚠ Some difference between C++ and R results\n")
    }
    
    # Performance comparison
    if (cpp_time > 0 && r_time > 0) {
      speedup <- r_time / cpp_time
      cat("  C++ speedup:", round(speedup, 2), "x\n")
    }
    
  }, error = function(e) {
    cat("✗ R implementation failed:", e$message, "\n")
  })
  
  # Test with different wavelets
  cat("\nTesting different wavelets...\n")
  wavelets <- c("la8", "haar")
  for (wf in wavelets) {
    tryCatch({
      result <- fast_wavelet_forward_1d(y, wf = wf, J = 3, boundary = "periodic")
      reconstructed <- fast_wavelet_inverse_1d(result$coeff, result$map)
      mse <- mean((y - reconstructed)^2)
      cat("  ", wf, ": MSE =", round(mse, 8), "\n")
    }, error = function(e) {
      cat("  ", wf, ": FAILED -", e$message, "\n")
    })
  }
  
  cat("\n=== C++ Wavelet Test Complete ===\n")
  return(TRUE)
}

# Run the test
test_wavelet_cpp()
