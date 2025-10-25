# Simple wavelet test to verify C++ implementation works
# This script tests the wavelet functions on a simple signal

# Load WICMAD package
library(devtools)
devtools::load_all()

# Test function
test_wavelets <- function() {
  cat("=== Simple Wavelet Test ===\n\n")
  
  # Create a simple test signal (dyadic length for R wavelets)
  set.seed(123)
  P <- 32  # 32 is 2^5, so dyadic
  t <- seq(0, 1, length.out = P)
  y <- sin(2*pi*t) + 0.1*rnorm(P)
  
  cat("Test signal created:\n")
  cat("- Length:", P, "\n")
  cat("- Range:", round(range(y), 3), "\n\n")
  
  # Test 1: Check if C++ functions are available
  cat("1. Checking C++ function availability...\n")
  if (exists("fast_wavelet_forward_1d")) {
    cat("✓ C++ function 'fast_wavelet_forward_1d' is available\n")
  } else {
    cat("✗ C++ function 'fast_wavelet_forward_1d' is NOT available\n")
    cat("  This means the C++ code hasn't been compiled yet.\n")
  }
  
  if (exists("fast_wavelet_inverse_1d")) {
    cat("✓ C++ function 'fast_wavelet_inverse_1d' is available\n")
  } else {
    cat("✗ C++ function 'fast_wavelet_inverse_1d' is NOT available\n")
  }
  
  # Test 2: Test R wavelet functions
  cat("\n2. Testing R wavelet functions...\n")
  start_time <- Sys.time()
  wt_result <- NULL
  tryCatch({
    wt_result <<- wt_forward_1d(y, wf = "la8", J = 3, boundary = "periodic")  # J=3 for 32=2^5
    r_time <- as.numeric(Sys.time() - start_time, units = "secs")
    cat("✓ R forward transform successful\n")
    cat("  Time:", round(r_time, 4), "seconds\n")
    cat("  Coefficients length:", length(wt_result$coeff), "\n")
  }, error = function(e) {
    cat("✗ R forward transform failed:", e$message, "\n")
    return(FALSE)
  })
  
  # Test inverse
  start_time <- Sys.time()
  tryCatch({
    reconstructed <- wt_inverse_1d(wt_result$coeff, wt_result$map)
    r_inv_time <- as.numeric(Sys.time() - start_time, units = "secs")
    cat("✓ R inverse transform successful\n")
    cat("  Time:", round(r_inv_time, 4), "seconds\n")
    
    # Check reconstruction accuracy
    mse <- mean((y - reconstructed)^2)
    cat("  Reconstruction MSE:", round(mse, 8), "\n")
    if (mse < 1e-10) {
      cat("  ✓ Perfect reconstruction!\n")
    } else if (mse < 1e-6) {
      cat("  ✓ Excellent reconstruction (numerical precision)\n")
    } else {
      cat("  ⚠ Some reconstruction error\n")
    }
  }, error = function(e) {
    cat("✗ R inverse transform failed:", e$message, "\n")
    return(FALSE)
  })
  
  # Test 3: Test C++ implementation if available
  if (exists("fast_wavelet_forward_1d")) {
    cat("\n3. Testing C++ wavelet functions...\n")
    start_time <- Sys.time()
    tryCatch({
      cpp_result <- fast_wavelet_forward_1d(y, wf = "la8", J = 3, boundary = "periodic")
      cpp_time <- as.numeric(Sys.time() - start_time, units = "secs")
      cat("✓ C++ forward transform successful\n")
      cat("  Time:", round(cpp_time, 4), "seconds\n")
      
      # Compare with R implementation if available
      if (!is.null(wt_result)) {
        coeff_diff <- mean(abs(wt_result$coeff - cpp_result$coeff))
        cat("  Coefficient difference from R:", round(coeff_diff, 8), "\n")
        
        if (coeff_diff < 1e-6) {
          cat("  ✓ C++ and R results match!\n")
        } else {
          cat("  ⚠ Some difference between C++ and R results\n")
        }
        
        # Performance comparison
        if (exists("r_time") && r_time > 0 && cpp_time > 0) {
          speedup <- r_time / cpp_time
          cat("  C++ speedup:", round(speedup, 2), "x\n")
        }
      } else {
        cat("  R implementation failed, cannot compare\n")
      }
      
    }, error = function(e) {
      cat("✗ C++ implementation failed:", e$message, "\n")
      cat("  Falling back to R implementation\n")
    })
  } else {
    cat("\n3. C++ functions not available - skipping C++ test\n")
  }
  
  # Test 4: Test with different wavelets
  cat("\n4. Testing different wavelets...\n")
  wavelets <- c("la8", "haar")
  for (wf in wavelets) {
    tryCatch({
      result <- wt_forward_1d(y, wf = wf, J = 3, boundary = "periodic")
      reconstructed <- wt_inverse_1d(result$coeff, result$map)
      mse <- mean((y - reconstructed)^2)
      cat("  ", wf, ": MSE =", round(mse, 8), "\n")
    }, error = function(e) {
      cat("  ", wf, ": FAILED -", e$message, "\n")
    })
  }
  
  cat("\n=== Test Complete ===\n")
  return(TRUE)
}

# Run the test
test_wavelets()
