# Julia scripts for running WICMAD experiments

This directory mirrors the dataset launchers from the original R codebase with lightweight Julia equivalents. The entry point is [`run_wicmad.jl`](run_wicmad.jl), which loads a UCR-style time-series dataset, constructs an imbalanced anomaly-detection split, and executes the Julia `wicmad` sampler. A dedicated ArrowHead runner lives in [`arrowhead/run_arrowhead.jl`](arrowhead/run_arrowhead.jl) for parity with the original `scripts/arrowhead` workflow.

## Usage

```bash
julia --project=julia julia/scripts/run_wicmad.jl \
  --dataset ArrowHead \
  --data-dir /path/to/UCRArchive2018/ArrowHead \
  --anomaly-ratio 0.10 \
  --reveal-ratio 0.15 \
  --n-iter 6000 \
  --burn 3000 \
  --thin 5
```

### Required arguments
- `--dataset` – dataset prefix used by the `.ts` files (for example `ArrowHead` expects `ArrowHead_TRAIN.ts` and `ArrowHead_TEST.ts`).
- `--data-dir` – directory containing the dataset files.

### Optional arguments
- `--anomaly-ratio` – target anomaly share when subsampling the majority/minority classes (default `0.10`).
- `--reveal-ratio` – fraction of anomalies to reveal (pin) during sampling (default `0.15`).
- `--n-iter`, `--burn`, `--thin` – MCMC configuration, matching the Julia `wicmad` keyword arguments.
- `--seed` – RNG seed used for subsampling, revealing anomalies, and the sampler itself.
- `--mean-intercept` – enable the optional intercept updates.
- `--no-diagnostics` – skip allocating diagnostic arrays to save memory.
- `--output` – if provided, writes a small TOML summary (ARI, clusters, anomaly fraction) to the given file.

The script activates the local Julia project, so dependencies declared in `julia/Project.toml` are available automatically.

## ArrowHead quick start

The Julia port of the ArrowHead workflow mirrors the R script located in `scripts/arrowhead`. It bakes in dataset-specific defaults (such as `n_iter = 8000` and `warmup = 500`) while still exposing flags for experimentation.

```bash
julia --project=julia julia/scripts/arrowhead/run_arrowhead.jl \
  --data-dir /path/to/UCRArchive2018/ArrowHead \
  --metrics results/arrowhead_metrics.toml
```

Key options:

- `--data-dir` – folder that contains `ArrowHead_TRAIN.ts` (and optionally `ArrowHead_TEST.ts`). Defaults to `julia/data/ArrowHead` if it exists.
- `--metrics` – optional TOML output summarising ARI, confusion matrix, and per-class precision/recall/F1.
- Sampling controls (`--n-iter`, `--burn`, `--thin`, `--warmup`, `--mean-intercept`, `--no-diagnostics`) align with the generic runner.
- Dataset preparation flags (`--anomaly-ratio`, `--reveal-ratio`, `--seed`) behave identically to `run_wicmad.jl` but target the majority/minority structure of ArrowHead specifically.

## Helper utilities

`common.jl` implements dataset loaders and sampling utilities:
- `load_ucr_dataset` – loads train/test `.ts` files and concatenates them.
- `prepare_anomaly_dataset` – converts the multiclass labels into a binary anomaly-detection setup with a configurable anomaly proportion.
- `summarize_dataset` and `default_time_index` – convenience helpers used by the main script.

These helpers are factored so additional dataset-specific launchers can reuse them.
