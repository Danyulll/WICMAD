#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(WICMAD)
})

# -----------------------------------------------------------------------------
# Global simulation controls
# -----------------------------------------------------------------------------
P_use <- 128L
iterations <- 10000L
burnin <- 2000L
thin <- 5L
use_parallel <- FALSE
warmup_iters <- 500L

output_root <- file.path("simstudy_runs", format(Sys.time(), "%Y%m%d_%H%M%S"))
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

message("Writing outputs to: ", normalizePath(output_root, winslash = "/", mustWork = FALSE))

model_configs <- list(
  list(
    key = "no_mean_intercept_bias_kernel",
    label = "No mean intercept + kernel bias",
    mean_intercept = FALSE,
    kernel_mode = "bias_only"
  ),
  list(
    key = "mean_intercept_bias_kernel",
    label = "Mean intercept + kernel bias",
    mean_intercept = TRUE,
    kernel_mode = "bias_only"
  ),
  list(
    key = "mean_intercept_no_bias",
    label = "Mean intercept + no kernel bias",
    mean_intercept = TRUE,
    kernel_mode = "no_bias"
  )
)

# -----------------------------------------------------------------------------
# Helper functions lifted from the simulation study vignette
# -----------------------------------------------------------------------------

dist_matrix <- function(x) outer(x, x, function(a, b) abs(a - b))

k_se <- function(x, ell, sf2) {
  D <- dist_matrix(x)
  sf2 * exp(-(D^2) / (2 * ell^2))
}

k_periodic <- function(x, per, ellp, sf2) {
  D <- dist_matrix(x)
  sf2 * exp(-2 * (sin(pi * D / per)^2) / (ellp^2))
}

k_quasiperiodic <- function(x, ell1 = 0.15, sf2_1 = 1.0, per = 0.30, ellp = 0.30,
                            sf2_per = 1.0, ell2 = 0.60, sf2_2 = 0.4) {
  k1 <- k_se(x, ell1, sf2_1) * k_periodic(x, per, ellp, sf2_per)
  k2 <- k_se(x, ell2, sf2_2)
  k1 + k2
}

f_base <- function(tt) 0.6 * sin(2 * pi * 2 * tt) + 0.25 * cos(2 * pi * 5 * tt) + 0.1 * tt
f_base_flat <- function(tt) 0.15 * sin(2 * pi * 2 * tt) + 0.08 * cos(2 * pi * 5 * tt) + 0.02 * tt

gp_sample <- function(grid, mean_fun, K, noise_sd = 0.05, jitter = 1e-6) {
  P <- length(grid)
  m <- mean_fun(grid)
  Sigma <- K + diag(noise_sd^2, P) + diag(jitter, P)
  L <- chol(Sigma)
  z <- rnorm(P)
  as.numeric(m + t(L) %*% z)
}

make_normal_curves_gp <- function(n = 40, P = 256, kernel_fun = k_quasiperiodic,
                                  noise_sd = 0.05, mean_fun = f_base) {
  grid <- seq(0, 1, length.out = P)
  K <- kernel_fun(grid)
  purrr::map_dfr(
    seq_len(n),
    function(i) tibble::tibble(id = i, t = grid, x = gp_sample(grid, mean_fun, K, noise_sd))
  )
}

make_multivariate_normals_gp <- function(N = 40, M = 3, P = 256, Q = 2,
                                         kernel_fun = k_quasiperiodic, noise_sd = 0.05,
                                         A = matrix(c(1.0, 0.7,
                                                      0.6, 1.0,
                                                      0.8, 0.4),
                                                    nrow = 3, byrow = TRUE)) {
  stopifnot(M == nrow(A), Q == ncol(A))
  grid <- seq(0, 1, length.out = P)
  K <- kernel_fun(grid)
  U <- replicate(Q, gp_sample(grid, f_base, K, noise_sd = 0.0))
  purrr::map_dfr(seq_len(N), function(i) {
    eps <- matrix(rnorm(P * M, 0, noise_sd), nrow = P, ncol = M)
    X <- U %*% t(A) + eps
    tibble::tibble(
      id = i,
      t = rep(grid, times = M),
      x = as.vector(X),
      channel = factor(rep(paste0("ch", seq_len(M)), each = P),
                       levels = paste0("ch", seq_len(M)))
    )
  })
}

gen_staerman_anomaly <- function(type, grid, scale_ref = 1.0, interval_frac = 0.10) {
  P <- length(grid)
  dx <- (max(grid) - min(grid)) / (P - 1)
  if (type == "isolated") {
    amp <- runif(1, 3, 4) * 2.5 * scale_ref
    sgn <- sample(c(-1, 1), 1)
    i0 <- sample(3:(P - 2), 1)
    w <- dx * runif(1, 0.3, 0.8)
    A <- amp * sgn
    A * exp(- (grid - grid[i0])^2 / (2 * w^2))
  } else if (type == "mag1") {
    rep(runif(1, 12, 15) * scale_ref * sample(c(-1, 1), 1), P)
  } else if (type == "mag2") {
    len <- max(2L, round(interval_frac * P))
    start <- sample(1:(P - len + 1), 1)
    lev <- runif(1, 10, 15) * scale_ref * sample(c(-1, 1), 1)
    v <- rep(0, P)
    idx <- start:(start + len - 1)
    s <- seq(0, 1, length.out = length(idx))
    ramp <- 0.5 * (1 - cos(pi * s))
    v[idx] <- lev * ramp
    v
  } else if (type == "shape") {
    u4 <- runif(1, 0.2, 2.0)
    3.0 * scale_ref * sin(2 * pi * u4 * grid)
  } else {
    stop("type must be one of: isolated, mag1, mag2, shape")
  }
}

num_deriv1 <- function(x, t) {
  P <- length(x)
  v <- numeric(P)
  v[1] <- (x[2] - x[1]) / (t[2] - t[1])
  v[P] <- (x[P] - x[P - 1]) / (t[P] - t[P - 1])
  if (P > 2) v[2:(P - 1)] <- (x[3:P] - x[1:(P - 2)]) / (t[3:P] - t[1:(P - 2)])
  v
}

num_deriv2 <- function(x, t) {
  P <- length(x)
  v <- numeric(P)
  dt <- mean(diff(t))
  if (P >= 3) v[2:(P - 1)] <- (x[3:P] - 2 * x[2:(P - 1)] + x[1:(P - 2)]) / (dt^2)
  v[1] <- v[2]
  v[P] <- v[P - 1]
  v
}

plot_functions_before <- function(df, title = "") {
  if (!"channel" %in% names(df)) df <- mutate(df, channel = factor("ch1", levels = "ch1"))
  ggplot(df, aes(t, x, group = interaction(id, channel), color = label)) +
    geom_line(alpha = 0.55, linewidth = 0.55) +
    scale_color_manual(values = c("Normal" = "#1f77b4", "Anomaly" = "#d62728")) +
    facet_wrap(~channel, ncol = 1, scales = "free_y") +
    labs(title = title, x = "t", y = "x(t)", color = "") +
    theme_minimal(base_size = 12)
}

with_kernel_bias_mode <- function(mode = c("default", "bias_only", "bias_enabled", "no_bias"), expr) {
  mode <- match.arg(mode)
  ns <- asNamespace("WICMAD")
  mk_orig <- get("make_kernels", envir = ns)
  restore <- function() {
    unlockBinding("make_kernels", ns)
    assign("make_kernels", mk_orig, envir = ns)
    lockBinding("make_kernels", ns)
  }
  new_fun <- switch(
    mode,
    bias_only = function(add_bias_variants = TRUE) {
      kernels <- mk_orig(TRUE)
      bias_kernels <- Filter(function(kc) isTRUE(grepl("\\+Bias$", kc$name)), kernels)
      if (length(bias_kernels) == 0) {
        stop("No bias-augmented kernels available in make_kernels().")
      }
      bias_kernels
    },
    bias_enabled = function(add_bias_variants = TRUE) mk_orig(TRUE),
    no_bias = function(add_bias_variants = TRUE) mk_orig(FALSE),
    default = mk_orig
  )
  unlockBinding("make_kernels", ns)
  assign("make_kernels", new_fun, envir = ns)
  lockBinding("make_kernels", ns)
  on.exit(restore(), add = TRUE)
  force(expr)
}

uni_to_model_raw <- function(df) {
  ids <- sort(unique(df$id))
  tvec <- df |> arrange(t) |> pull(t) |> unique()
  Y <- lapply(ids, function(i) df |> filter(id == i) |> arrange(t) |> pull(x) |> matrix(ncol = 1L))
  by_id <- df |> distinct(id, label) |> arrange(id)
  list(Y = Y, t = tvec, ids = ids, is_anom = ifelse(by_id$label == "Anomaly", 1L, 0L))
}

uni_to_model_deriv12 <- function(df, scale_each = FALSE) {
  ids <- sort(unique(df$id))
  tvec <- df |> arrange(t) |> pull(t) |> unique()
  Y <- lapply(ids, function(i) {
    xi <- df |> filter(id == i) |> arrange(t) |> pull(x)
    d1 <- num_deriv1(xi, tvec)
    d2 <- num_deriv2(xi, tvec)
    if (scale_each) {
      s1 <- sd(d1)
      if (is.finite(s1) && s1 > 0) d1 <- d1 / s1
      s2 <- sd(d2)
      if (is.finite(s2) && s2 > 0) d2 <- d2 / s2
    }
    cbind(xi, d1, d2)
  })
  by_id <- df |> distinct(id, label) |> arrange(id)
  list(Y = Y, t = tvec, ids = ids, is_anom = ifelse(by_id$label == "Anomaly", 1L, 0L))
}

svd_components_dataset <- function(df, r = 4) {
  ids <- df |> distinct(id) |> arrange(id) |> pull(id)
  tvec <- df |> arrange(t) |> pull(t) |> unique()
  N <- length(ids)
  P <- length(tvec)
  X <- matrix(NA_real_, nrow = N, ncol = P)
  lab_by_id <- df |> distinct(id, label) |> arrange(id) |> pull(label)
  for (ii in seq_along(ids)) {
    xi <- df |> filter(id == ids[ii]) |> arrange(t) |> pull(x)
    stopifnot(length(xi) == P)
    X[ii, ] <- xi
  }
  Xc <- scale(X, center = TRUE, scale = FALSE)
  sv <- svd(Xc)
  r <- min(r, length(sv$d))
  U <- sv$u[, 1:r, drop = FALSE]
  D <- diag(sv$d[1:r], r, r)
  V <- sv$v[, 1:r, drop = FALSE]
  scores <- U %*% D
  comps_list <- vector("list", N)
  for (i in seq_len(N)) {
    comp_mat <- matrix(NA_real_, nrow = P, ncol = r)
    for (k in seq_len(r)) comp_mat[, k] <- scores[i, k] * V[, k]
    comps_list[[i]] <- tibble::tibble(
      id = ids[i],
      t = rep(tvec, times = r),
      x = as.vector(comp_mat),
      channel = factor(rep(paste0("comp", seq_len(r)), each = P),
                       levels = paste0("comp", seq_len(r))),
      label = factor(rep(lab_by_id[i], each = P * r), levels = c("Normal", "Anomaly"))
    )
  }
  df_comps <- bind_rows(comps_list)
  Y <- lapply(seq_len(N), function(i) {
    df_comps |> filter(id == ids[i]) |> arrange(t) |>
      select(t, channel, x) |>
      pivot_wider(names_from = channel, values_from = x) |>
      select(all_of(paste0("comp", seq_len(r)))) |>
      as.matrix()
  })
  list(df = df_comps, Y = Y, t = tvec, ids = ids, labels = lab_by_id)
}

plot_clustered_inline <- function(Y, t, z_hat, title = "Clustered curves (Dahl)") {
  stopifnot(length(Y) == length(z_hat))
  df_list <- lapply(seq_along(Y), function(i) {
    mat <- as.matrix(Y[[i]])
    P_i <- nrow(mat)
    M_i <- ncol(mat)
    tibble::tibble(
      id = i,
      t = rep(t, times = M_i),
      value = as.vector(mat),
      channel = factor(rep(paste0("ch", seq_len(M_i)), each = P_i)),
      cluster = factor(z_hat[i])
    )
  })
  dfl <- bind_rows(df_list)
  cols <- scales::hue_pal()(length(unique(dfl$cluster)))
  ggplot(dfl, aes(t, value, group = interaction(id, channel), color = cluster)) +
    geom_line(alpha = 0.7, linewidth = 0.5) +
    scale_color_manual(values = cols) +
    facet_wrap(~channel, ncol = 1, scales = "free_y") +
    labs(x = "t", y = "x(t)", title = title, color = "cluster") +
    theme_minimal(base_size = 12)
}

# Normal/anomaly binning utilities ------------------------------------------------

determine_normal_anomaly_clusters <- function(
    z_hat,
    reveal_idx = integer(0),
    min_size = 5L,
    min_prop = NULL) {
  z_hat <- as.integer(z_hat)
  N <- length(z_hat)
  stopifnot(N >= 1L, all(z_hat >= 1L))
  cl_ids <- sort(unique(z_hat))
  cl_sizes <- vapply(cl_ids, function(k) sum(z_hat == k), integer(1))
  names(cl_sizes) <- as.character(cl_ids)
  thr_abs <- as.integer(min_size)
  thr_prop <- if (!is.null(min_prop)) ceiling(min_prop * N) else 0L
  thr <- max(thr_abs, thr_prop)
  is_small <- cl_sizes <= thr
  eligible <- cl_ids[!is_small]
  reveal_idx <- unique(as.integer(reveal_idx))
  reveal_idx <- reveal_idx[reveal_idx >= 1 & reveal_idx <= N]
  has_revealed <- length(reveal_idx) > 0
  choose_normal <- function(cands) {
    if (length(cands) == 0L) return(integer(0))
    if (has_revealed) {
      rev_counts <- vapply(cands, function(k) sum(z_hat[reveal_idx] == k), integer(1))
      names(rev_counts) <- as.character(cands)
      tied <- cands[rev_counts == max(rev_counts)]
      if (length(tied) > 1L) tied[which.max(cl_sizes[as.character(tied)])] else tied
    } else {
      cands[which.max(cl_sizes[as.character(cands)])]
    }
  }
  normal_cluster_id <- choose_normal(eligible)
  if (length(normal_cluster_id) == 0L) normal_cluster_id <- choose_normal(cl_ids)
  anomaly_cluster_ids <- setdiff(cl_ids, normal_cluster_id)
  pred_anom <- as.integer(z_hat %in% anomaly_cluster_ids)
  list(
    normal_cluster_id = as.integer(normal_cluster_id),
    anomaly_cluster_ids = as.integer(anomaly_cluster_ids),
    pred_anom = pred_anom,
    cluster_sizes = cl_sizes
  )
}

confusion_after <- function(df_long, z_hat, reveal_idx = integer(0),
                            min_size = 5L, min_prop = NULL,
                            dataset_name = NULL, make_plot = TRUE) {
  bin <- determine_normal_anomaly_clusters(
    z_hat = z_hat,
    reveal_idx = reveal_idx,
    min_size = min_size,
    min_prop = min_prop
  )
  truth_lab <- df_long |>
    dplyr::distinct(id, label) |>
    dplyr::arrange(id) |>
    dplyr::pull(label)
  true_anom <- as.integer(truth_lab == "Anomaly")
  stopifnot(length(true_anom) == length(bin$pred_anom))
  confmat <- table(
    Truth = factor(ifelse(true_anom == 1, "Anomaly", "Normal"),
                   levels = c("Normal", "Anomaly")),
    Pred = factor(ifelse(bin$pred_anom == 1, "Anomaly", "Normal"),
                  levels = c("Normal", "Anomaly"))
  )
  TP <- confmat["Anomaly", "Anomaly"]
  TN <- confmat["Normal", "Normal"]
  FP <- confmat["Normal", "Anomaly"]
  FN <- confmat["Anomaly", "Normal"]
  acc <- (TP + TN) / sum(confmat)
  prec <- if ((TP + FP) > 0) TP / (TP + FP) else NA
  rec <- if ((TP + FN) > 0) TP / (TP + FN) else NA
  f1 <- if (!is.na(prec) && !is.na(rec) && (prec + rec) > 0) 2 * prec * rec / (prec + rec) else NA
  ttl <- if (is.null(dataset_name)) "Confusion Matrix" else paste("Confusion Matrix —", dataset_name)
  p <- NULL
  if (isTRUE(make_plot)) {
    cm_df <- as.data.frame(confmat)
    p <- ggplot(cm_df, aes(x = Pred, y = Truth, fill = Freq)) +
      geom_tile(color = "white") +
      geom_text(aes(label = Freq), fontface = "bold") +
      scale_fill_gradient(low = "#f7fbff", high = "#08306b") +
      labs(title = ttl, x = "Predicted", y = "True") +
      theme_minimal()
  }
  list(confmat = confmat,
       metrics = list(acc = acc, prec = prec, rec = rec, f1 = f1),
       bin = bin,
       plot = p)
}

# -----------------------------------------------------------------------------
# Data generators (adapted for P_use = 128)
# -----------------------------------------------------------------------------

generate_isolated <- function(P = P_use, N_curves = 50L) {
  set.seed(82123)
  grid <- seq(0, 1, length.out = P)
  K <- k_quasiperiodic(grid)
  df_iso <- purrr::map_dfr(seq_len(N_curves), function(i) {
    xi <- gp_sample(grid, f_base_flat, K, noise_sd = 0.02)
    tibble::tibble(id = i, t = grid, x = xi)
  })
  n_anom <- max(1L, floor(0.05 * N_curves))
  ids_anom <- sort(sample(seq_len(N_curves), n_anom, replace = FALSE))
  for (i in ids_anom) {
    i0 <- sample(3:(P - 2), 1)
    w <- (1 / (P - 1)) * runif(1, 0.3, 0.8)
    A <- runif(1, 8, 12) * sample(c(-1, 1), 1)
    bump <- A * exp(-(grid - grid[i0])^2 / (2 * w^2))
    df_iso$x[df_iso$id == i] <- df_iso$x[df_iso$id == i] + bump
  }
  df_iso |>
    mutate(label = if_else(id %in% ids_anom, "Anomaly", "Normal"),
           label = factor(label, levels = c("Normal", "Anomaly")))
}

generate_mag1 <- function(P = P_use, N_curves = 40L) {
  set.seed(8746)
  grid <- seq(0, 1, length.out = P)
  df_mag1 <- make_normal_curves_gp(n = N_curves, P = P)
  n_anom <- max(1L, floor(0.05 * N_curves))
  ids_anom <- sort(sample(seq_len(N_curves), n_anom, replace = FALSE))
  for (i in ids_anom) {
    df_mag1$x[df_mag1$id == i] <- df_mag1$x[df_mag1$id == i] + gen_staerman_anomaly("mag1", grid, 1.0)
  }
  df_mag1 |>
    mutate(label = if_else(id %in% ids_anom, "Anomaly", "Normal"),
           label = factor(label, levels = c("Normal", "Anomaly")))
}

generate_mag2 <- function(P = P_use, N_curves = 60L) {
  set.seed(473423)
  grid <- seq(0, 1, length.out = P)
  df_mag2 <- make_normal_curves_gp(n = N_curves, P = P)
  n_anom <- max(1L, floor(0.05 * N_curves))
  ids_anom <- sort(sample(seq_len(N_curves), n_anom, replace = FALSE))
  for (i in ids_anom) {
    df_mag2$x[df_mag2$id == i] <- df_mag2$x[df_mag2$id == i] +
      gen_staerman_anomaly("mag2", grid, 1.0, interval_frac = 0.10)
  }
  df_mag2 |>
    mutate(label = if_else(id %in% ids_anom, "Anomaly", "Normal"),
           label = factor(label, levels = c("Normal", "Anomaly")))
}

generate_shape <- function(P = P_use, N_curves = 60L) {
  set.seed(9623)
  grid <- seq(0, 1, length.out = P)
  df_shape <- make_normal_curves_gp(n = N_curves, P = P)
  n_anom <- max(1L, floor(0.05 * N_curves))
  ids_anom <- sort(sample(seq_len(N_curves), n_anom, replace = FALSE))
  for (i in ids_anom) {
    df_shape$x[df_shape$id == i] <- df_shape$x[df_shape$id == i] + gen_staerman_anomaly("shape", grid, 1.0)
  }
  df_shape |>
    mutate(label = if_else(id %in% ids_anom, "Anomaly", "Normal"),
           label = factor(label, levels = c("Normal", "Anomaly")))
}

generate_mv_one <- function(P = P_use, N = 40L, M = 3L) {
  set.seed(9307)
  base_mv <- make_multivariate_normals_gp(N = N, M = M, P = P)
  ids_pool <- sort(unique(base_mv$id))
  n_anom <- max(1L, floor(0.05 * length(ids_pool)))
  n_anom <- min(n_anom, length(ids_pool))
  ids_anom <- sort(sample(ids_pool, size = n_anom, replace = FALSE))
  plan <- tidyr::expand_grid(id = ids_pool, channel = levels(base_mv$channel)) |>
    mutate(type = NA_character_)
  for (i in ids_anom) {
    ch <- sample(levels(base_mv$channel), 1)
    aty <- sample(c("isolated", "mag1", "mag2", "shape"), 1)
    plan$type[plan$id == i & plan$channel == ch] <- aty
  }
  grid <- base_mv |> arrange(t) |> pull(t) |> unique()
  df_mv_one <- base_mv |>
    left_join(plan, by = c("id", "channel")) |>
    group_by(id, channel) |>
    group_modify(~{
      if (!is.na(.x$type[1])) .x$x <- .x$x + gen_staerman_anomaly(.x$type[1], grid, 1.0)
      .x
    }) |>
    ungroup() |>
    mutate(label = if_else(id %in% ids_anom, "Anomaly", "Normal"),
           label = factor(label, levels = c("Normal", "Anomaly")))
  df_mv_one
}

generate_mv_two <- function(P = P_use, N = 40L, M = 3L) {
  set.seed(21983)
  base_mv <- make_multivariate_normals_gp(N = N, M = M, P = P)
  ids_pool <- sort(unique(base_mv$id))
  n_anom <- max(1L, floor(0.05 * length(ids_pool)))
  ids_anom <- sort(sample(ids_pool, size = n_anom, replace = FALSE))
  plan <- tidyr::expand_grid(id = ids_pool, channel = levels(base_mv$channel)) |>
    mutate(type = NA_character_)
  for (i in ids_anom) {
    chs <- sample(levels(base_mv$channel), 2, replace = FALSE)
    aty <- sample(c("isolated", "mag1", "mag2", "shape"), 1)
    plan$type[plan$id == i & plan$channel %in% chs] <- aty
  }
  grid <- base_mv |> arrange(t) |> pull(t) |> unique()
  base_mv |>
    left_join(plan, by = c("id", "channel")) |>
    group_by(id, channel) |>
    group_modify(~{
      if (!is.na(.x$type[1])) .x$x <- .x$x + gen_staerman_anomaly(.x$type[1], grid, 1.0)
      .x
    }) |>
    ungroup() |>
    mutate(label = if_else(id %in% ids_anom, "Anomaly", "Normal"),
           label = factor(label, levels = c("Normal", "Anomaly")))
}

generate_mv_three <- function(P = P_use, N = 40L, M = 3L) {
  set.seed(48571)
  base_mv <- make_multivariate_normals_gp(N = N, M = M, P = P)
  ids_pool <- sort(unique(base_mv$id))
  n_anom <- max(1L, floor(0.05 * length(ids_pool)))
  ids_anom <- sort(sample(ids_pool, size = n_anom, replace = FALSE))
  plan <- tidyr::expand_grid(id = ids_pool, channel = levels(base_mv$channel)) |>
    mutate(type = NA_character_)
  for (i in ids_anom) {
    atys <- sample(c("isolated", "mag1", "mag2", "shape"), 3, replace = FALSE)
    plan$type[plan$id == i & plan$channel == "ch1"] <- atys[1]
    plan$type[plan$id == i & plan$channel == "ch2"] <- atys[2]
    plan$type[plan$id == i & plan$channel == "ch3"] <- atys[3]
  }
  grid <- base_mv |> arrange(t) |> pull(t) |> unique()
  base_mv |>
    left_join(plan, by = c("id", "channel")) |>
    group_by(id, channel) |>
    group_modify(~{
      if (!is.na(.x$type[1])) .x$x <- .x$x + gen_staerman_anomaly(.x$type[1], grid, 1.0)
      .x
    }) |>
    ungroup() |>
    mutate(label = if_else(id %in% ids_anom, "Anomaly", "Normal"),
           label = factor(label, levels = c("Normal", "Anomaly")))
}

# -----------------------------------------------------------------------------
# Task constructors
# -----------------------------------------------------------------------------

sample_revealed_idx <- function(df, prop = 0.15) {
  normal_ids <- df |> distinct(id, label) |> filter(label == "Normal") |> pull(id)
  if (length(normal_ids) == 0) return(integer(0))
  n_pick <- max(1L, round(prop * length(normal_ids)))
  sort(sample(normal_ids, size = n_pick, replace = FALSE))
}

create_univariate_tasks <- function(df, base_id, base_title, K_init = 8L,
                                    min_size = 3L, min_prop = 0.05) {
  revealed_idx <- sample_revealed_idx(df)
  ds_raw <- uni_to_model_raw(df)
  svd_ds <- svd_components_dataset(df, r = 4)
  df_deriv <- df |>
    group_by(id) |>
    arrange(t, .by_group = TRUE) |>
    mutate(dx = num_deriv1(x, t), ddx = num_deriv2(x, t)) |>
    ungroup() |>
    pivot_longer(cols = c(x, dx, ddx), names_to = "channel", values_to = "val") |>
    mutate(channel = factor(channel, levels = c("x", "dx", "ddx")),
           x = val) |>
    select(-val)
  ds_deriv <- uni_to_model_deriv12(df, scale_each = FALSE)
  list(
    list(
      id = paste0(base_id, "_raw"),
      display_name = paste0(base_title, " — Raw"),
      plot_df = df,
      Y = ds_raw$Y,
      t = ds_raw$t,
      df_truth = df,
      revealed_idx = revealed_idx,
      K_init = K_init,
      min_size = min_size,
      min_prop = min_prop
    ),
    list(
      id = paste0(base_id, "_svd4"),
      display_name = paste0(base_title, " — SVD(4)"),
      plot_df = svd_ds$df,
      Y = svd_ds$Y,
      t = svd_ds$t,
      df_truth = df,
      revealed_idx = revealed_idx,
      K_init = K_init,
      min_size = min_size,
      min_prop = min_prop
    ),
    list(
      id = paste0(base_id, "_deriv12"),
      display_name = paste0(base_title, " — Deriv(1,2)"),
      plot_df = df_deriv,
      Y = ds_deriv$Y,
      t = ds_deriv$t,
      df_truth = df,
      revealed_idx = revealed_idx,
      K_init = K_init,
      min_size = min_size,
      min_prop = min_prop
    )
  )
}

create_mv_task <- function(df, id, title, K_init = 5L, min_size = 2L, min_prop = 0.05) {
  ids <- sort(unique(df$id))
  chs <- sort(unique(df$channel))
  tvec <- df |> arrange(t) |> pull(t) |> unique()
  Y <- lapply(ids, function(i) {
    df |> filter(id == i) |> arrange(t) |>
      select(t, channel, x) |>
      pivot_wider(names_from = channel, values_from = x) |>
      select(all_of(chs)) |>
      as.matrix()
  })
  list(
    id = id,
    display_name = title,
    plot_df = df,
    Y = Y,
    t = tvec,
    df_truth = df,
    revealed_idx = sample_revealed_idx(df),
    K_init = K_init,
    min_size = min_size,
    min_prop = min_prop
  )
}

# -----------------------------------------------------------------------------
# Assemble task list
# -----------------------------------------------------------------------------

all_tasks <- list()

all_tasks <- append(all_tasks, create_univariate_tasks(generate_isolated(),
                                                     base_id = "model1_isolated",
                                                     base_title = "Model 1 — Isolated",
                                                     K_init = 10L,
                                                     min_size = 3L,
                                                     min_prop = 0.05))

all_tasks <- append(all_tasks, create_univariate_tasks(generate_mag1(),
                                                     base_id = "model2_mag1",
                                                     base_title = "Model 2 — Magnitude I",
                                                     K_init = 8L,
                                                     min_size = 3L,
                                                     min_prop = NULL))

all_tasks <- append(all_tasks, create_univariate_tasks(generate_mag2(),
                                                     base_id = "model3_mag2",
                                                     base_title = "Model 3 — Magnitude II",
                                                     K_init = 8L,
                                                     min_size = 3L,
                                                     min_prop = 0.05))

all_tasks <- append(all_tasks, create_univariate_tasks(generate_shape(),
                                                     base_id = "model4_shape",
                                                     base_title = "Model 4 — Shape",
                                                     K_init = 8L,
                                                     min_size = 3L,
                                                     min_prop = 0.05))

all_tasks <- append(all_tasks,
  list(create_mv_task(generate_mv_one(), id = "mv_one_channel",
                      title = "MV: One Channel", K_init = 5L,
                      min_size = 2L, min_prop = 0.05)))

all_tasks <- append(all_tasks,
  list(create_mv_task(generate_mv_two(), id = "mv_two_channels",
                      title = "MV: Two Channels (Same)", K_init = 5L,
                      min_size = 2L, min_prop = 0.05)))

all_tasks <- append(all_tasks,
  list(create_mv_task(generate_mv_three(), id = "mv_three_channels",
                      title = "MV: Three Mixed", K_init = 5L,
                      min_size = 2L, min_prop = 0.05)))

dataset_count <- length(all_tasks)
all_jobs <- list()

for (i in seq_along(all_tasks)) {
  dataset <- all_tasks[[i]]
  dataset_dir <- file.path(output_root, sprintf("%02d_%s", i, dataset$id))
  for (cfg in model_configs) {
    job <- list(
      dataset_id = dataset$id,
      dataset_label = dataset$display_name,
      plot_df = dataset$plot_df,
      Y = dataset$Y,
      t = dataset$t,
      df_truth = dataset$df_truth,
      revealed_idx = dataset$revealed_idx,
      K_init = dataset$K_init,
      min_size = dataset$min_size,
      min_prop = dataset$min_prop,
      dataset_dir = dataset_dir,
      output_dir = file.path(dataset_dir, cfg$key),
      config_key = cfg$key,
      config_label = cfg$label,
      mean_intercept = cfg$mean_intercept,
      kernel_mode = cfg$kernel_mode,
      job_id = paste(dataset$id, cfg$key, sep = "__")
    )
    job$display_name <- paste0(dataset$display_name, " — ", cfg$label)
    all_jobs <- append(all_jobs, list(job))
  }
}

total_jobs <- length(all_jobs)
threads_needed <- total_jobs
message(sprintf(
  "Total datasets: %d. Model configurations per dataset: %d.",
  dataset_count, length(model_configs)
))

# -----------------------------------------------------------------------------
# Runner
# -----------------------------------------------------------------------------

run_dataset_task <- function(task) {
  dir.create(task$dataset_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(task$output_dir, recursive = TRUE, showWarnings = FALSE)
  before_plot <- plot_functions_before(task$plot_df,
                                       paste0(task$display_name, " — dataset to be clustered"))
  before_path <- file.path(task$output_dir, "01_before.png")
  ggplot2::ggsave(before_path, before_plot, dpi = 150, width = 9, height = 6, units = "in")

  message(sprintf(
    "[%s] Starting wicmad() run (n_iter=%d, burn=%d, thin=%d, K_init=%d, mean_intercept=%s, kernel_mode=%s)",
    task$job_id, iterations, burnin, thin, task$K_init,
    ifelse(task$mean_intercept, "TRUE", "FALSE"), task$kernel_mode
  ))

  run_model <- function() {
    wicmad(
      Y = task$Y,
      t = task$t,
      mean_intercept = task$mean_intercept,
      n_iter = iterations,
      burn = burnin,
      thin = thin,
      K_init = task$K_init,
      revealed_idx = task$revealed_idx,
      warmup_iters = warmup_iters,
      unpin = FALSE,
      diagnostics = TRUE
    )
  }

  res <- switch(
    task$kernel_mode,
    bias_only = with_kernel_bias_mode("bias_only", run_model()),
    bias_enabled = with_kernel_bias_mode("bias_enabled", run_model()),
    no_bias = with_kernel_bias_mode("no_bias", run_model()),
    with_kernel_bias_mode("default", run_model())
  )

  res_path <- file.path(task$output_dir, "wicmad_result.rds")
  saveRDS(res, res_path)

  dahl <- dahl_from_res(res)
  clustered_plot <- plot_clustered_inline(task$Y, task$t, dahl$z_hat,
                                          paste0(task$display_name, " — clustered (Dahl)"))
  cluster_path <- file.path(task$output_dir, "02_clustered.png")
  ggplot2::ggsave(cluster_path, clustered_plot, dpi = 150, width = 9, height = 6, units = "in")

  cm <- confusion_after(task$df_truth, dahl$z_hat,
                        reveal_idx = task$revealed_idx,
                        min_size = task$min_size,
                        min_prop = task$min_prop,
                        dataset_name = task$display_name,
                        make_plot = TRUE)
  if (!is.null(cm$plot)) {
    conf_path <- file.path(task$output_dir, "03_confusion.png")
    ggplot2::ggsave(conf_path, cm$plot, dpi = 150, width = 6, height = 5, units = "in")
  }

  diag_dir <- file.path(task$output_dir, "diagnostics")
  save_diagnostics_from_res(res, Y = task$Y, t = task$t, out_dir = diag_dir,
                            dpi = 150, width = 9, height = 6, units = "in")

  list(
    dataset_id = task$dataset_id,
    dataset_label = task$dataset_label,
    config_key = task$config_key,
    config_label = task$config_label,
    mean_intercept = task$mean_intercept,
    kernel_mode = task$kernel_mode,
    job_id = task$job_id,
    output_dir = task$output_dir,
    result_path = res_path,
    confusion = cm$metrics
  )
}

threads_to_use <- min(threads_needed, parallel::detectCores(logical = TRUE))
if (threads_to_use < threads_needed) {
  warning(sprintf("Only %d cores detected; datasets will be processed on %d workers instead of %d.",
                  parallel::detectCores(logical = TRUE), threads_to_use, threads_needed))
}

if (threads_to_use <= 1L) {
  results <- lapply(all_jobs, run_dataset_task)
} else {
  cl <- parallel::makeCluster(threads_to_use)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterEvalQ(cl, {
    suppressPackageStartupMessages({
      library(dplyr)
      library(tidyr)
      library(purrr)
      library(ggplot2)
      library(WICMAD)
    })
  })
  parallel::clusterExport(cl,
    varlist = c("run_dataset_task", "iterations", "burnin", "thin", "warmup_iters",
                "plot_functions_before", "plot_clustered_inline", "confusion_after",
                "determine_normal_anomaly_clusters", "output_root", "use_parallel",
                "with_kernel_bias_mode"),
    envir = environment()
  )
  results <- parallel::parLapply(cl, all_jobs, run_dataset_task)
}

message("All dataset runs completed.")
summary_path <- file.path(output_root, "summary_metrics.csv")
summary_df <- purrr::map_dfr(results, function(x) {
  tibble::tibble(
    dataset_id = x$dataset_id,
    dataset_label = x$dataset_label,
    config_key = x$config_key,
    config_label = x$config_label,
    mean_intercept = x$mean_intercept,
    kernel_mode = x$kernel_mode,
    job_id = x$job_id,
    output_dir = x$output_dir,
    accuracy = x$confusion$acc,
    precision = x$confusion$prec,
    recall = x$confusion$rec,
    f1 = x$confusion$f1
  )
})
utils::write.csv(summary_df, summary_path, row.names = FALSE)
message("Saved summary metrics to ", normalizePath(summary_path, winslash = "/", mustWork = FALSE))
