# ArrowHead Subsample Test Summary

## Test Scripts Created

### 1. `test_arrowhead_subsample.R` - Comprehensive Test
- **Purpose**: Test WICMAD algorithm on ArrowHead subsample
- **Features**:
  - Loads 20 samples from ArrowHead dataset
  - Tests wavelet functions for correctness
  - Runs WICMAD clustering (50 iterations)
  - Measures performance metrics
  - Checks C++ implementation availability

### 2. `simple_wavelet_test.R` - Basic Wavelet Test
- **Purpose**: Simple test of wavelet functions
- **Features**:
  - Tests on synthetic signal
  - Compares R vs C++ implementations
  - Tests different wavelets (la8, haar)
  - Measures reconstruction accuracy

## What the Tests Verify

### ✅ Wavelet Function Correctness
1. **Forward Transform**: `wt_forward_1d()` works correctly
2. **Inverse Transform**: `wt_inverse_1d()` reconstructs perfectly
3. **Multiple Wavelets**: Tests la8 and haar wavelets
4. **Reconstruction Accuracy**: MSE < 1e-6 for perfect reconstruction

### ✅ C++ Implementation
1. **Function Availability**: Checks if C++ functions are compiled
2. **Correctness**: Compares C++ vs R results
3. **Performance**: Measures speedup achieved
4. **Fallback**: Ensures R implementation works if C++ fails

### ✅ WICMAD Integration
1. **Data Loading**: ArrowHead subsample loads correctly
2. **Preprocessing**: Interpolation to p=32 works
3. **Clustering**: WICMAD runs successfully
4. **Parameters**: alpha_prior=c(5,1) produces reasonable clusters

## Expected Results

### If C++ is Compiled:
```
✓ C++ function 'fast_wavelet_forward_1d' is available
✓ C++ function 'fast_wavelet_inverse_1d' is available
✓ C++ forward transform successful
✓ C++ and R results match!
C++ speedup: 3-5x
```

### If C++ is NOT Compiled:
```
✗ C++ function 'fast_wavelet_forward_1d' is NOT available
✓ R forward transform successful
✓ R inverse transform successful
✓ Perfect reconstruction!
```

## Performance Expectations

### Current Performance (R only):
- Wavelet operations: ~1.78s self time (20% of total)
- Total analysis: ~25.7s for 100 iterations

### Expected Performance (C++):
- Wavelet operations: ~0.45s self time (3-5x speedup)
- Total analysis: ~24.4s for 100 iterations (5.2% improvement)

## Test Execution

### Method 1: Direct R Execution
```bash
R --vanilla --slave -e "source('scripts/arrowhead/simple_wavelet_test.R')"
```

### Method 2: Batch File
```bash
.\run_subsample_test.bat
```

### Method 3: R Interactive
```r
source('scripts/arrowhead/simple_wavelet_test.R')
```

## Troubleshooting

### If C++ Functions Not Available:
1. **Compile C++ code**: Run `Rcpp::compileAttributes()`
2. **Rebuild package**: `devtools::load_all()`
3. **Check compilation**: Look for errors in compilation

### If Wavelet Functions Fail:
1. **Check dependencies**: Ensure `waveslim` package is installed
2. **Check data**: Verify ArrowHead data files exist
3. **Check parameters**: Ensure p=32 is reasonable

### If WICMAD Fails:
1. **Check data format**: Ensure series are matrices
2. **Check parameters**: Verify alpha_prior is reasonable
3. **Check iterations**: Start with minimal iterations

## Success Criteria

### ✅ Test Passes If:
- Wavelet functions work correctly
- Reconstruction MSE < 1e-6
- WICMAD completes successfully
- Performance metrics are collected

### ❌ Test Fails If:
- Wavelet functions throw errors
- Reconstruction is inaccurate
- WICMAD algorithm fails
- C++ compilation errors

## Next Steps After Testing

1. **If C++ works**: Run full profiler to measure speedup
2. **If C++ fails**: Debug compilation issues
3. **If R works**: Verify fallback implementation
4. **If all fails**: Check package dependencies

## Files to Check

- `R/RcppExports.R` - Should contain C++ function exports
- `src/fast_icm_loglik_curve.cpp` - C++ implementation
- `R/wavelets.R` - R wrapper functions
- `scripts/arrowhead/` - Test scripts

The test scripts are ready to verify that the wavelet C++ implementation works correctly on real ArrowHead data!
