# WICMAD.jl

A Julia port of the WICMAD (Wavelet-ICM Anomaly Detector) sampler. The implementation
mirrors the original R package's Dirichlet Process mixture sampler with Besov-shrunk
wavelet blocks and intrinsic coregionalization kernels.

## Getting started

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()

using WICMAD
```

## Running the sampler

```julia
Y = [randn(64, 2) for _ in 1:10]
t = range(0, 1; length=64)
result = WICMAD.wicmad(Y, collect(t); n_iter=100, burn=50, thin=5, diagnostics=false)
```

The sampler returns a named tuple containing cluster allocations, kernel selections,
posterior draws for the concentration parameter, and diagnostic summaries.

## Notes

* Wavelet transforms rely on the `Wavelets.jl` package; ensure the observation length is a power of two.
* Diagnostics are intentionally lightweight compared with the R package but track
  the key global quantities (occupied clusters, concentration parameter, log-likelihood).
