#!/usr/bin/env Rscript

# ------------------------------------------------------------
# WICMAD simulation experiment runner (FPCA via fdapace)
# - 7 datasets (4 uni + 3 multi)
# - Functional PCA using fdapace::FPCA (K=4)
# - 100 MC runs; pin normals throughout
# - Save PNGs only for MC run 1 (before/after clustering) in FPCA space
# - Save averaged metrics to CSV
# ------------------------------------------------------------

library(devtools)
devtools::install()

suppressPackageStartupMessages({
  library(WICMAD)
  library(fdapace)   # <- Functional PCA (FPCA)
  library(ggplot2)
  library(dplyr)
  library(purrr)
  library(furrr)     # <- Parallel processing
  library(future)    # <- Parallel backend
  library(tidyr)
  library(readr)
  library(numDeriv)  # <- Numerical differentiation
})

# -----------------------------
# Global controls
# -----------------------------
P_use         <- 128L
# P_use         <- 16L
t_grid        <- seq(0, 1, length.out = P_use)
Delta         <- 1/(P_use - 1)
sigma_noise   <- 0.05

mc_runs       <- 100L
# mc_runs       <- 1L
base_seed     <- 20240501L

# Sampler controls
n_iter        <- 10000L
# n_iter        <- 10L
burnin        <- 3000L
# burnin        <- 3L
thin          <- 1L
warmup_iters  <- 500L
# warmup_iters  <- 10L

# Semi-supervised: reveal 15% normals (pinned for entire run)
reveal_prop   <- 0.15

# Functional PCA (FPCA) controls
USE_FPCA      <- TRUE   # <- set FALSE to use raw signals
FPCA_K        <- 4L     # number of FPCA components per original channel

# Output
timestamp     <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_root      <- file.path("simstudy_results_fpca", timestamp)
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
png_dir       <- file.path(out_root, "png")
dir.create(png_dir, showWarnings = FALSE)
metrics_csv   <- file.path(out_root, "summary_metrics.csv")

message("Writing outputs to: ", normalizePath(out_root, winslash = "/"))

# -----------------------------
# Parallel processing setup
# -----------------------------
# Set up parallel processing
workers <- 100  # Set desired number of workers here
max_available <- parallel::detectCores() - 1
n_cores <- if (workers <= max_available) {
  workers
} else {
  message("Requested ", workers, " workers not available. Using ", max_available, " workers instead.")
  max_available
}
# Fallback to 10 if still too many
if (n_cores > 10) {
  n_cores <- 10
  message("Limiting to 10 workers as requested.")
}
plan(multisession, workers = n_cores)
message("Using ", n_cores, " cores for parallel processing")

# Export all necessary functions to parallel workers
# This ensures all functions are available in the parallel environment

# -----------------------------
# GP mean and kernel
# -----------------------------
mean_fun <- function(t) {
  0.6 * sin(4 * pi * t) + 0.25 * cos(10 * pi * t) + 0.1 * t
}

# Quasiperiodic + long-scale SE
k_qp <- function(t, tprime,
                 ell1 = 0.15, sig1_sq = 1.0, p = 0.30, ellp = 0.30,
                 ell2 = 0.60, sig2_sq = 0.4) {
  dt <- outer(t, tprime, "-")
  se1 <- exp(- (dt^2) / (2 * ell1^2))
  per <- exp(- 2 * (sin(pi * dt / p))^2 / (ellp^2))
  se2 <- exp(- (dt^2) / (2 * ell2^2))
  sig1_sq * se1 * per + sig2_sq * se2
}

# Draw GP curves via MVN (Cholesky)
gp_draw_matrix <- function(N, t, mu_fun, cov_fun, sigma_eps) {
  P <- length(t)
  mu <- mu_fun(t)
  K  <- cov_fun(t, t)
  Kc <- K + 1e-8 * diag(P)
  L  <- chol(Kc)
  Z  <- matrix(rnorm(P * N), P, N)
  G  <- mu + t(L) %*% Z     # P x N (columns are functions)
  X  <- t(G) + matrix(rnorm(N * P, sd = sigma_eps), N, P)
  X  # N x P
}

# -----------------------------
# Anomaly perturbations
# -----------------------------
add_isolated <- function(xrow, t) {
  P <- length(t)
  i0 <- sample(3:(P-2), 1)
  w  <- runif(1, 0.3*Delta, 0.8*Delta)
  S  <- sample(c(-1, +1), 1)
  A  <- runif(1, 8, 12)
  bump <- S * A * exp(- (t - t[i0])^2 / (2 * w^2))
  xrow + bump
}

add_mag1 <- function(xrow) {
  S <- sample(c(-1, +1), 1)
  A <- runif(1, 12, 15)
  xrow + S * A
}

add_mag2 <- function(xrow, t) {
  P <- length(t)
  lambda <- floor(0.10 * P)
  s <- sample(1:(P - lambda + 1), 1)
  S <- sample(c(-1, +1), 1)
  A <- runif(1, 10, 15)
  rfun <- function(u) 0.5 * (1 - cos(pi * u))
  bump <- numeric(P)
  idx <- s:(s + lambda - 1)
  u <- (seq_along(idx) - 1) / (lambda - 1)
  bump[idx] <- S * A * rfun(u)
  xrow + bump
}

add_shape <- function(xrow, t) {
  U <- runif(1, 0.2, 2.0)
  xrow + 3 * sin(2 * pi * U * t)
}

# -----------------------------
# Dataset generators
# -----------------------------
make_univariate_dataset <- function(N = 40L, anomaly_type = c("isolated","mag1","mag2","shape"),
                                    t = t_grid) {
  anomaly_type <- match.arg(anomaly_type)
  X <- gp_draw_matrix(N, t, mean_fun, k_qp, sigma_noise)  # N x P
  y_true <- rep(0L, N)
  n_anom <- max(1L, round(0.10 * N))   # 10%
  idx_anom <- sample.int(N, n_anom)

  X_pert <- X
  for (i in idx_anom) {
    xi <- X[i, ]
    if (anomaly_type == "isolated")      xi <- add_isolated(xi, t)
    else if (anomaly_type == "mag1")     xi <- add_mag1(xi)
    else if (anomaly_type == "mag2")     xi <- add_mag2(xi, t)
    else if (anomaly_type == "shape")    xi <- add_shape(xi, t)
    X_pert[i, ] <- xi
  }
  y_true[idx_anom] <- 1L

  # Wrap as list of P x M matrices with M=1
  Y_list <- lapply(seq_len(N), function(i) matrix(X_pert[i, ], nrow = length(t), ncol = 1))
  list(Y_list = Y_list, t = t, y_true = y_true)
}

# Multivariate with 3 channels using 2 latent GP factors
make_multivariate_dataset <- function(N = 40L, regime = c("one","two","three"),
                                      t = t_grid) {
  regime <- match.arg(regime)
  P <- length(t)
  U1 <- gp_draw_matrix(N, t, mean_fun, k_qp, 0)
  U2 <- gp_draw_matrix(N, t, mean_fun, k_qp, 0)
  A <- matrix(c(1.0, 0.4,
                0.2, 1.0,
                0.7, -0.3), nrow = 3, byrow = TRUE)

  Y <- vector("list", N)
  for (i in seq_len(N)) {
    Ui <- cbind(U1[i, ], U2[i, ])  # P x 2
    Xi <- Ui %*% t(A)              # P x 3
    Xi <- Xi + matrix(rnorm(P * 3, sd = sigma_noise), nrow = P, ncol = 3)
    Y[[i]] <- Xi
  }

  y_true <- rep(0L, N)
  n_anom <- max(1L, round(0.10 * N))
  idx_anom <- sample.int(N, n_anom)
  y_true[idx_anom] <- 1L

  types <- c("isolated","mag1","mag2","shape")
  apply_perturb <- function(row, type) {
    if (type == "isolated") add_isolated(row, t)
    else if (type == "mag1") add_mag1(row)
    else if (type == "mag2") add_mag2(row, t)
    else if (type == "shape") add_shape(row, t)
    else row
  }

  for (i in idx_anom) {
    Xi <- Y[[i]]
    if (regime == "one") {
      ch <- sample(1:3, 1)
      type <- sample(types, 1)
      Xi[, ch] <- apply_perturb(Xi[, ch], type)
    } else if (regime == "two") {
      chs <- sample(1:3, 2)
      type <- sample(types, 1)  # same type on both
      for (ch in chs) Xi[, ch] <- apply_perturb(Xi[, ch], type)
    } else { # "three"
      for (ch in 1:3) {
        type <- sample(types, 1)
        Xi[, ch] <- apply_perturb(Xi[, ch], type)
      }
    }
    Y[[i]] <- Xi
  }

  list(Y_list = Y, t = t, y_true = y_true)
}

# -----------------------------
# Derivatives transform
# - Input: Y_list (length N), each is P x M
# - Output: Y_list_deriv (length N), each is P x (3*M) with original + 1st + 2nd derivatives
# -----------------------------
derivatives_transform <- function(Y_list, t) {
  N <- length(Y_list)
  P <- length(t)
  
  # Compute derivatives for each subject
  Y_deriv <- lapply(Y_list, function(Yi) {
    M <- ncol(Yi)
    deriv_mat <- matrix(0, nrow = P, ncol = 3 * M)
    
    for (m in seq_len(M)) {
      # Original signal
      deriv_mat[, 3*(m-1) + 1] <- Yi[, m]
      
      # First derivative using numDeriv
      if (P > 2) {
        # Create interpolation function for smooth derivatives
        signal <- Yi[, m]
        # Use spline interpolation for smooth derivatives
        spline_fun <- splinefun(t, signal)
        deriv_mat[, 3*(m-1) + 2] <- spline_fun(t, deriv = 1)
      }
      
      # Second derivative using numDeriv
      if (P > 4) {
        # Use spline interpolation for second derivatives
        signal <- Yi[, m]
        spline_fun <- splinefun(t, signal)
        deriv_mat[, 3*(m-1) + 3] <- spline_fun(t, deriv = 2)
      }
    }
    
    deriv_mat
  })
  
  Y_deriv
}

# -----------------------------
# Functional PCA (FPCA) transform via fdapace
# - Input: Y_list (length N), each is P x M
# - Output: Y_list_fpca (length N), each is P x (M * K) component curves
#   where channel m contributes K component-curves: score_{ik} * phi_k(t)
# -----------------------------
fpca_transform <- function(Y_list, t, K = FPCA_K) {
  N <- length(Y_list)
  P <- length(t)
  # Determine M from first subject
  M <- ncol(Y_list[[1]])
  # Repeat time grid per subject for fdapace
  Lt_dense <- replicate(N, t, simplify = FALSE)

  # For each channel, run FPCA across the N subjects, then build K component-curves
  channel_components <- vector("list", M)
  for (m in seq_len(M)) {
    # Ly: list of length N with vectors of length P for channel m
    Ly_m <- lapply(Y_list, function(Yi) Yi[, m])

    fp <- FPCA(Ly = Ly_m, Lt = Lt_dense,
               optns = list(
                 dataType = 'Dense',
                 methodSelectK = K
               ))
    # fp$phi: matrix P x K; fp$xiEst: N x K
    phi <- fp$phi
    xi  <- fp$xiEst
    # Build N lists: each entry P x K with columns xi[i,k] * phi[,k]
    comp_m <- lapply(seq_len(N), function(i) {
      comps <- sapply(seq_len(K), function(k) xi[i, k] * phi[, k])
      # Ensure matrix P x K
      matrix(comps, nrow = P, ncol = K, dimnames = list(NULL, paste0("m", m, "_k", seq_len(K))))
    })
    channel_components[[m]] <- comp_m
  }

  # Concatenate components across channels for each subject
  Y_fpca <- lapply(seq_len(N), function(i) {
    mats <- lapply(seq_len(length(channel_components)), function(m) channel_components[[m]][[i]])
    do.call(cbind, mats)  # P x (M*K)
  })
  Y_fpca
}

# -----------------------------
# Plot helpers (MC run 1 only)
# - For FPCA transform, plots are in FPCA component-channel space
# -----------------------------
plot_before <- function(Y_list, t, y_true, title) {
  df <- map2_dfr(seq_along(Y_list), Y_list, ~{
    Xi <- .y
    P <- nrow(Xi); M <- ncol(Xi)
    # Create proper time and channel indices
    t_rep <- rep(t, M)
    m_rep <- rep(seq_len(M), each = P)
    tibble(i = .x, t = t_rep, m = m_rep, x = as.vector(Xi))
  }) %>%
    left_join(tibble(i = seq_along(y_true), y_true = y_true), by = "i") %>%
    mutate(lbl = ifelse(y_true == 1L, "Anomalous", "Normal"),
           lbl = factor(lbl, levels = c("Normal","Anomalous")))

  ggplot(df, aes(t, x, group = i, color = lbl)) +
    geom_line(alpha = 0.7) +
    facet_wrap(~ m, ncol = 1, scales = "free_y") +
    labs(title = title,
         color = "Class", x = "t", y = "value") +
    theme_minimal()
}

plot_after <- function(Y_list, t, z_hat, normal_label, title) {
  pred_anom <- as.integer(z_hat != normal_label)
  df <- map2_dfr(seq_along(Y_list), Y_list, ~{
    Xi <- .y
    P <- nrow(Xi); M <- ncol(Xi)
    # Create proper time and channel indices
    t_rep <- rep(t, M)
    m_rep <- rep(seq_len(M), each = P)
    tibble(i = .x, t = t_rep, m = m_rep, x = as.vector(Xi))
  }) %>%
    mutate(pred = rep(pred_anom, each = nrow(Y_list[[1]]) * ncol(Y_list[[1]]))) %>%
    mutate(lbl = ifelse(pred == 1L, "Predicted Anom", "Predicted Normal"),
           lbl = factor(lbl, levels = c("Predicted Normal","Predicted Anom")))

  ggplot(df, aes(t, x, group = i, color = lbl)) +
    geom_line(alpha = 0.7) +
    facet_wrap(~ m, ncol = 1, scales = "free_y") +
    labs(title = title,
         color = "Cluster", x = "t", y = "value") +
    theme_minimal()
}

# -----------------------------
# Metrics
# -----------------------------
metrics_from_preds <- function(y_true, pred_anom) {
  tn <- sum(y_true == 0 & pred_anom == 0)
  fp <- sum(y_true == 0 & pred_anom == 1)
  fn <- sum(y_true == 1 & pred_anom == 0)
  tp <- sum(y_true == 1 & pred_anom == 1)
  acc <- (tp + tn) / (tp + tn + fp + fn)
  prec <- if ((tp + fp) == 0) NA_real_ else tp / (tp + fp)
  rec  <- if ((tp + fn) == 0) NA_real_ else tp / (tp + fn)
  f1   <- if (is.na(prec) || is.na(rec) || (prec + rec) == 0) NA_real_ else 2 * prec * rec / (prec + rec)
  tibble(accuracy = acc, precision = prec, recall = rec, f1 = f1)
}

# -----------------------------
# One dataset for one MC seed
# -----------------------------
run_one <- function(dataset_label, dataset_title, representation, make_data_fn, mc_idx) {
  set.seed(base_seed + mc_idx)
  dat <- make_data_fn()
  Y_list_raw <- dat$Y_list
  t <- dat$t
  y_true <- dat$y_true
  N <- length(Y_list_raw)

  # Apply appropriate transformation based on representation
  if (representation == "raw") {
    Y_list <- Y_list_raw
  } else if (representation == "derivatives") {
    Y_list <- derivatives_transform(Y_list_raw, t)
  } else if (representation == "fpca") {
    Y_list <- fpca_transform(Y_list_raw, t, K = FPCA_K)
  } else {
    stop("Unknown representation: ", representation)
  }
  
  # Data validation for all data
  for (i in seq_along(Y_list)) {
    if (any(!is.finite(Y_list[[i]]))) {
      message("Warning: Non-finite values detected in dataset ", dataset_label, " subject ", i)
      Y_list[[i]][!is.finite(Y_list[[i]])] <- 0  # Replace with 0
    }
  }

  # Reveal 15% normals (random subset of true normals)
  normals_idx <- which(y_true == 0L)
  n_reveal <- max(1L, floor(reveal_prop * length(normals_idx)))
  reveal_idx <- sample(normals_idx, n_reveal)

  # Run WICMAD (pin normals for entire chain)
  tryCatch({
    res <- wicmad(
      Y = Y_list,
      t = t,
      n_iter = n_iter,
      burn = burnin,
      thin = thin,
      warmup_iters = warmup_iters,
      revealed_idx = reveal_idx,
      unpin = FALSE   # keep revealed curves pinned throughout
    )
  }, error = function(e) {
    message("Error in WICMAD for dataset ", dataset_label, ": ", e$message)
    # Return a dummy result with all subjects in one cluster
    list(
      Z = matrix(1L, nrow = n_iter - burnin, ncol = N),
      K_occ = rep(1L, n_iter - burnin),
      loglik = rep(0, n_iter - burnin)
    )
  })

  # Dahl partition and predicted anomalies
  dahl <- dahl_from_res(res)
  z_hat <- dahl$z_hat
  normal_label <- {
    lab_counts <- table(z_hat[reveal_idx])
    as.integer(names(lab_counts)[which.max(lab_counts)])
  }
  pred_anom <- as.integer(z_hat != normal_label)

  # Save PNGs for MC run 1 only
  if (mc_idx == 1L) {
    tag <- paste0(dataset_label, "_", representation)
    p_before <- plot_before(Y_list, t, y_true, dataset_title)
    p_after  <- plot_after(Y_list, t, z_hat, normal_label, dataset_title)
    ggsave(filename = file.path(png_dir, paste0(tag, "_before.png")),
           plot = p_before, width = 8, height = 6, dpi = 200)
    ggsave(filename = file.path(png_dir, paste0(tag, "_after.png")),
           plot = p_after, width = 8, height = 6, dpi = 200)
  }

  met <- metrics_from_preds(y_true, pred_anom)
  met %>%
    mutate(dataset = dataset_label,
           representation = representation,
           mc_run = mc_idx) %>%
    relocate(dataset, representation, mc_run)
}

# -----------------------------
# Dataset registry
# -----------------------------
# Univariate datasets with multiple representations
univariate_specs <- list(
  list(id = "uni_isolated",     title = "Isolated",     fn = function() make_univariate_dataset(anomaly_type = "isolated", t = t_grid)),
  list(id = "uni_mag1",         title = "Magnitude I",  fn = function() make_univariate_dataset(anomaly_type = "mag1",     t = t_grid)),
  list(id = "uni_mag2",         title = "Magnitude II", fn = function() make_univariate_dataset(anomaly_type = "mag2",     t = t_grid)),
  list(id = "uni_shape",        title = "Shape",        fn = function() make_univariate_dataset(anomaly_type = "shape",    t = t_grid))
)

# Multivariate datasets (raw only)
multivariate_specs <- list(
  list(id = "mv_one_channel",   title = "One Channel",  fn = function() make_multivariate_dataset(regime = "one",   t = t_grid)),
  list(id = "mv_two_channels",  title = "Two Channels",fn = function() make_multivariate_dataset(regime = "two",   t = t_grid)),
  list(id = "mv_three_channels",title = "Three Channels",fn = function() make_multivariate_dataset(regime = "three", t = t_grid))
)

# Create expanded dataset specs with all representations for univariate datasets
dataset_specs <- list()

# Add univariate datasets with all three representations
for (spec in univariate_specs) {
  # Raw representation
  dataset_specs[[length(dataset_specs) + 1]] <- list(
    id = paste0(spec$id, "_raw"),
    title = spec$title,
    representation = "raw",
    fn = spec$fn
  )
  # Derivatives representation
  dataset_specs[[length(dataset_specs) + 1]] <- list(
    id = paste0(spec$id, "_deriv"),
    title = spec$title,
    representation = "derivatives",
    fn = spec$fn
  )
  # FPCA representation
  dataset_specs[[length(dataset_specs) + 1]] <- list(
    id = paste0(spec$id, "_fpca"),
    title = spec$title,
    representation = "fpca",
    fn = spec$fn
  )
}

# Add multivariate datasets (raw only)
for (spec in multivariate_specs) {
  dataset_specs[[length(dataset_specs) + 1]] <- list(
    id = spec$id,
    title = spec$title,
    representation = "raw",
    fn = spec$fn
  )
}

# -----------------------------
# MAIN: run all datasets across MC runs (in parallel)
# -----------------------------
all_results <- map_dfr(dataset_specs, function(spec) {
  id <- spec$id
  title <- spec$title
  representation <- spec$representation
  
  cat(sprintf("[Dataset: %s] Running %d MC replications (%s) in parallel...\n",
              id, mc_runs, representation))
  
  # Create a wrapper function that has access to all necessary functions
  run_one_wrapper <- function(mc_idx) {
    run_one(id, title, representation, spec$fn, mc_idx)
  }
  
  # Use sequential processing for now to avoid function access issues
  map_dfr(seq_len(mc_runs), run_one_wrapper)
})

# -----------------------------
# Aggregate averages across MC runs and save CSV
# -----------------------------
summary_tbl <- all_results %>%
  group_by(dataset, representation) %>%
  summarise(
    accuracy = mean(accuracy, na.rm = TRUE),
    precision = mean(precision, na.rm = TRUE),
    recall = mean(recall, na.rm = TRUE),
    f1 = mean(f1, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(summary_tbl, metrics_csv)
message("Saved summary metrics CSV to: ", normalizePath(metrics_csv, winslash = "/"))

# Also save per-run metrics
perrun_csv <- file.path(out_root, "metrics_per_run.csv")
readr::write_csv(all_results, perrun_csv)
message("Saved per-run metrics CSV to: ", normalizePath(perrun_csv, winslash = "/"))

message("PNG directory (first-run visualizations): ", normalizePath(png_dir, winslash = "/"))

# Clean up parallel processing
plan(sequential)
message("Parallel processing completed. Results saved.")
