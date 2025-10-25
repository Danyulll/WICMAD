#!/usr/bin/env Rscript

# Compilation script for C++ wavelet implementation
# This script compiles the C++ code and updates the RcppExports

cat("=== Compiling C++ Wavelet Implementation ===\n\n")

# Load required libraries
library(Rcpp)
library(RcppEigen)

# Set working directory to project root
setwd("../../")

cat("Current working directory:", getwd(), "\n")

# Compile attributes to generate RcppExports
cat("1. Compiling Rcpp attributes...\n")
tryCatch({
  Rcpp::compileAttributes()
  cat("✓ Rcpp attributes compiled successfully\n")
}, error = function(e) {
  cat("✗ Failed to compile Rcpp attributes:", e$message, "\n")
  stop("Compilation failed")
})

# Check if RcppExports.R was updated
if (file.exists("R/RcppExports.R")) {
  cat("✓ RcppExports.R updated\n")
} else {
  cat("✗ RcppExports.R not found\n")
}

# Check if the new functions are in RcppExports
if (file.exists("R/RcppExports.R")) {
  exports_content <- readLines("R/RcppExports.R")
  if (any(grepl("fast_wavelet_forward_1d", exports_content))) {
    cat("✓ fast_wavelet_forward_1d found in RcppExports\n")
  } else {
    cat("✗ fast_wavelet_forward_1d NOT found in RcppExports\n")
  }
  
  if (any(grepl("fast_wavelet_inverse_1d", exports_content))) {
    cat("✓ fast_wavelet_inverse_1d found in RcppExports\n")
  } else {
    cat("✗ fast_wavelet_inverse_1d NOT found in RcppExports\n")
  }
}

# Try to load the package to test compilation
cat("\n2. Testing package loading...\n")
tryCatch({
  devtools::load_all()
  cat("✓ Package loaded successfully\n")
  
  # Test if functions are available
  if (exists("fast_wavelet_forward_1d")) {
    cat("✓ fast_wavelet_forward_1d is available\n")
  } else {
    cat("✗ fast_wavelet_forward_1d is NOT available\n")
  }
  
  if (exists("fast_wavelet_inverse_1d")) {
    cat("✓ fast_wavelet_inverse_1d is available\n")
  } else {
    cat("✗ fast_wavelet_inverse_1d is NOT available\n")
  }
  
}, error = function(e) {
  cat("✗ Package loading failed:", e$message, "\n")
  cat("This might be due to compilation errors in the C++ code.\n")
})

cat("\n=== Compilation Complete ===\n")
cat("Next steps:\n")
cat("1. Run the test script: Rscript test_wavelet_cpp.R\n")
cat("2. If successful, run the profiling analysis to measure speedup\n")
