#' WICMAD: Wavelet-ICM DP-GP Sampler with Diagnostics
#'
#' Main MCMC sampler for the WICMAD (Wavelet-ICM Anomaly Detector) model, implementing a Dirichlet Process
#' mixture model for multivariate functional data with wavelet shrinkage and intrinsic
#' coregionalization modeling.
#'
#' @param Y List of length N, where each element is a P×M matrix representing a multivariate
#'   functional curve with P time points and M channels
#' @param t Numeric vector of length P or P×d matrix of time coordinates
#' @param n_iter Integer number of MCMC iterations (default: 6000)
#' @param burn Integer number of burn-in iterations (default: 3000)
#' @param thin Integer thinning interval for saved samples (default: 5)
#' @param alpha_prior Numeric vector of length 2 giving shape and rate parameters for
#'   the Gamma prior on the DP concentration parameter α (default: c(10, 1))
#' @param wf Character string specifying wavelet family (default: "la8")
#' @param J Integer number of wavelet decomposition levels or NULL for automatic
#' @param boundary Character string specifying boundary condition (default: "periodic")
#' @param mh_step_L Numeric proposal standard deviation for L matrix updates (default: 0.03)
#' @param mh_step_eta Numeric proposal standard deviation for η updates (default: 0.10)
#' @param mh_step_tauB Numeric proposal standard deviation for τ_B updates (default: 0.15)
#' @param revealed_idx Integer vector of curve indices to pin to cluster 1 during warmup
#' @param K_init Integer initial number of clusters (default: 5)
#' @param warmup_iters Integer number of warmup iterations for revealed curves (default: 100)
#' @param unpin Logical; if TRUE, release revealed curves after warmup (default: FALSE)
#' @param kappa_pi Numeric sparsity parameter for Besov prior (default: 0.6)
#' @param c2 Numeric decay parameter for level-specific sparsity (default: 1.0)
#' @param tau_pi Numeric concentration parameter for Beta prior (default: 40)
#' @param a_sig Numeric shape parameter for σ²_m prior (default: 2.5)
#' @param b_sig Numeric rate parameter for σ²_m prior (default: 0.02)
#' @param a_tau Numeric shape parameter for τ_σ prior (default: 2.0)
#' @param b_tau Numeric rate parameter for τ_σ prior (default: 2.0)
#' @param a_eta Numeric shape parameter for η prior (default: 2)
#' @param b_eta Numeric rate parameter for η prior (default: 0.1)
#' @param diagnostics Logical; if TRUE, collect diagnostic information (default: TRUE)
#' @param track_ids Integer vector of curve indices to track in diagnostics
#' @param monitor_levels Character vector of wavelet levels to monitor in diagnostics
#' @param mean_intercept Logical; if TRUE, include per-channel intercept in each cluster
#'   while disabling bias-augmented kernel variants (default: FALSE)
#'
#' @return List containing:
#'   \item{Z}{Integer matrix of cluster assignments (S×N)}
#'   \item{K_occ}{Integer vector of number of occupied clusters per iteration}
#'   \item{loglik}{Numeric vector of log-likelihood values}
#'   \item{diagnostics}{List of diagnostic information (if diagnostics=TRUE)}
#'   \item{params}{List of cluster-specific parameters}
#'   \item{meta}{List of metadata including dimensions and settings}
#'
#' @details The WICMAD model combines:
#' \itemize{
#'   \item \strong{Dirichlet Process clustering} for flexible number of clusters
#'   \item \strong{Wavelet shrinkage} using Besov priors for sparsity
#'   \item \strong{Intrinsic Coregionalization Model} for multivariate dependencies
#'   \item \strong{Gaussian Process} kernels for temporal/spatial correlation
#' }
#'
#' The MCMC algorithm includes:
#' \itemize{
#'   \item Cluster assignment updates using slice sampling for the Dirichlet Process
#'   \item Wavelet coefficient updates with spike-and-slab priors
#'   \item Kernel parameter updates via Metropolis-Hastings
#'   \item ICM parameter updates for multivariate structure
#' }
#'
#' @examples
#' # Basic usage with synthetic data
#' set.seed(123)
#' N <- 20; P <- 64; M <- 2
#' Y <- lapply(1:N, function(i) {
#'   t <- seq(0, 1, length.out = P)
#'   matrix(c(sin(2*pi*t) + rnorm(P, 0, 0.1),
#'            cos(2*pi*t) + rnorm(P, 0, 0.1)), ncol = M)
#' })
#' t <- seq(0, 1, length.out = P)
#' 
#' # Run MCMC
#' result <- wicmad(Y, t, n_iter = 1000, burn = 500, thin = 2)
#' 
#' # Extract cluster assignments
#' cluster_assignments <- result$Z
#' 
#' @export
wicmad <- function(
    Y, t,
    n_iter = 6000, burn = 3000, thin = 5,
    alpha_prior = c(10, 1),
    wf = "la8", J = NULL, boundary = "periodic",
    mh_step_L = 0.03, mh_step_eta = 0.10, mh_step_tauB = 0.15,
    revealed_idx = integer(0),
    K_init = 5,
    warmup_iters = 100,
    unpin = FALSE,
    kappa_pi = 0.6, c2 = 1.0, tau_pi = 40,
    a_sig = 2.5, b_sig = 0.02, a_tau = 2.0, b_tau = 2.0,
    a_eta = 2, b_eta = 0.1,
    diagnostics = TRUE,
    track_ids = NULL,
    monitor_levels = NULL,
    mean_intercept = FALSE
) {
  # ===== PACKAGE LOADING =====
  suppressPackageStartupMessages({
    requireNamespace("waveslim")
    requireNamespace("mvtnorm")
    requireNamespace("MASS")
    requireNamespace("invgamma")
  })
  
  # ===== PARAMETER VALIDATION AND SETUP =====
  # Guards so we always have post-burn iterations and (if feasible) a saved draw
  if (!is.numeric(n_iter) || n_iter < 2) stop("n_iter must be >= 2.")
  if (!is.numeric(thin) || thin < 1) thin <- 1L
  n_iter <- as.integer(n_iter); burn <- as.integer(burn); thin <- as.integer(thin)
  if (burn >= n_iter) {
    warning(sprintf("burn (%d) >= n_iter (%d). Resetting burn to n_iter - max(2, thin).", burn, n_iter))
    burn <- max(0L, n_iter - max(2L, thin))
  }
  if ((n_iter - burn) < thin) {
    new_thin <- max(1L, n_iter - burn)
    warning(sprintf("thin (%d) > n_iter - burn (%d). Using thin=%d so at least one draw is saved.",
                    thin, n_iter - burn, new_thin))
    thin <- new_thin
  }
  
  # ===== DATA DIMENSIONS AND TIME NORMALIZATION =====
  N <- length(Y); P <- nrow(Y[[1]]); M <- ncol(Y[[1]]); J <- ensure_dyadic_J(P, J)
  t <- normalize_t(t, P); t <- scale_t01(t)
  
  # ===== REVEALED CURVES HANDLING =====
  if (length(revealed_idx) > 0) {
    if (!unpin) {
      message(sprintf("Pinning %d revealed curves to cluster 1 for ALL %d iterations (unpin=FALSE).",
                      length(revealed_idx), n_iter))
    } else if (warmup_iters > 0) {
      message(sprintf("Warm-up: pinning %d revealed curves to cluster 1 for %d iterations (then release).",
                      length(revealed_idx), warmup_iters))
    } else {
      message(sprintf("Note: warmup_iters=0 and unpin=TRUE, so %d revealed curves are NOT pinned.",
                      length(revealed_idx)))
    }
  }
  
  # ===== INITIALIZATION: KERNELS, DP PARAMETERS, AND CLUSTERS =====
  kernels <- make_kernels(add_bias_variants = !isTRUE(mean_intercept))
  alpha <- rgamma(1, shape = alpha_prior[1], rate = alpha_prior[2])
  v     <- rbeta(K_init, 1, alpha)
  pi    <- stick_to_pi(v)
  K     <- length(v)
  
  # ===== CLUSTER PARAMETER INITIALIZATION =====
  params <- vector("list", K)
  for (k in 1:K) {
    params[[k]] <- draw_new_cluster_params(M, P, t, kernels, wf, J, boundary)
    params <- ensure_complete_cache(params, kernels, k, t, M)
  }
  
  # ===== INITIAL CLUSTER ASSIGNMENTS =====
  z <- sample.int(K, N, replace = TRUE)
  if (length(revealed_idx)) z[revealed_idx] <- 1L
  
  # ===== OUTPUT STORAGE INITIALIZATION =====
  keep <- max(0L, floor((n_iter - burn) / thin))
  Z_s <- if (keep > 0) matrix(NA_integer_, keep, N) else matrix(NA_integer_, 0, N)
  alpha_s <- if (keep > 0) numeric(keep) else numeric(0)
  kern_s  <- if (keep > 0) integer(keep) else integer(0)
  K_s     <- if (keep > 0) integer(keep) else integer(0)
  loglik_s<- if (keep > 0) numeric(keep) else numeric(0)
  
  # ===== HELPER FUNCTIONS =====
  # Adjusted Rand Index - helper function for computing ARI between clusterings
  adj_rand_index <- function(z1, z2) {
    tab <- table(z1, z2)
    sum_comb <- sum(choose(tab, 2))
    a <- rowSums(tab); b <- colSums(tab)
    sum_a <- sum(choose(a, 2)); sum_b <- sum(choose(b, 2))
    tot <- choose(length(z1), 2)
    (sum_comb - (sum_a*sum_b)/tot) / (0.5*(sum_a+sum_b) - (sum_a*sum_b)/tot)
  }

  # ===== DIAGNOSTICS INITIALIZATION =====
  diag <- NULL
  if (isTRUE(diagnostics)) {
    zeros <- matrix(0, P, M)
    tmp_fw <- wt_forward_mat(zeros, wf, J, boundary)
    lev_all <- names(tmp_fw[[1]]$map$idx)
    det_lev <- lev_all[grepl("^d", lev_all)]
    if (is.null(monitor_levels)) {
      j_low  <- det_lev[1]
      j_mid  <- det_lev[max(1, round(length(det_lev)/2))]
      j_high <- det_lev[length(det_lev)]
      monitor_levels <- unique(c(j_low, j_mid, j_high))
    }
    if (is.null(track_ids)) track_ids <- seq_len(min(5, N))
    
    diag <- list(
      meta = list(P=P, M=M, J=J, monitor_levels=monitor_levels, track_ids=track_ids),
      global = list(K_occ=if (keep>0) rep(NA_real_, keep) else numeric(0),
                    alpha= if (keep>0) rep(NA_real_, keep) else numeric(0),
                    loglik=if (keep>0) rep(NA_real_, keep) else numeric(0),
                    K_occ_all=if (n_iter >= 1) rep(NA_real_, n_iter) else numeric(0)),
      clusters = list(A=list(), B=list()),
      ari = if (keep>1) rep(NA_real_, keep-1) else numeric(0),
      ari_all = if (n_iter >= 1) rep(NA_real_, n_iter) else numeric(0)
    )
  }
  
  sidx <- 0L
  Y <- lapply(Y, as_num_mat)
  precomp_all <- precompute_wavelets(Y, wf, J, boundary)
  message(sprintf("Starting MCMC (N=%d, P=%d, M=%d, iters=%d; burn=%d, thin=%d, keep=%d)",
                  N, P, M, n_iter, burn, thin, keep))
  
  # ===== BIAS AND MEAN FUNCTION HELPER FUNCTIONS =====
  # Add bias terms to mean function - helper function
  add_bias_to_mu <- function(mu_wave, bias_vec) {
    if (!isTRUE(mean_intercept)) return(mu_wave)
    if (length(bias_vec) == 0) return(mu_wave)
    mu_wave + matrix(rep(bias_vec, each = nrow(mu_wave)), nrow(mu_wave), length(bias_vec))
  }

  # Convert bias terms to wavelet coefficients - helper function
  bias_to_wavelet <- function(bias_vec) {
    if (!isTRUE(mean_intercept)) return(NULL)
    if (length(bias_vec) == 0L) return(NULL)
    if (!any(abs(bias_vec) > 1e-12)) return(NULL)
    bias_mat <- matrix(rep(bias_vec, each = P), nrow = P, ncol = M)
    fw <- wt_forward_mat(bias_mat, wf, J, boundary)
    lapply(fw, function(comp) comp$coeff)
  }

  # Compute log-likelihood for curve in cluster - helper function
  ll_curve_k <- function(k, Yi, mu_k) {
    kc <- kernels[[ params[[k]]$kern_idx ]]
    kp <- params[[k]]$thetas[[ params[[k]]$kern_idx ]]
    cache <- .build_icm_cache(t, kc, kp, params[[k]]$L, params[[k]]$eta, params[[k]]$tau_B, params[[k]]$cache)
    params[[k]]$cache <<- cache
    y_resid <- Yi - mu_k
    fast_icm_loglik_curve(y_resid, cache)
  }

  # Sample channel-specific bias terms - helper function
  sample_bias <- function(k, idx, mu_wave) {
    if (!isTRUE(mean_intercept) || length(idx) == 0L) return(params[[k]]$bias)
    bias_cur <- params[[k]]$bias
    bias_sd <- 5
    prior_prec <- if (bias_sd > 0) 1 / (bias_sd^2) else 0
    for (m in seq_len(M)) {
      eta_m <- params[[k]]$eta[m]
      if (!is.finite(eta_m) || eta_m <= 0) eta_m <- 1e-6
      resid_sum <- 0
      for (ii in idx) {
        resid_sum <- resid_sum + sum(Y[[ii]][, m] - mu_wave[, m])
      }
      like_prec <- length(idx) * P / eta_m
      post_prec <- like_prec + prior_prec
      if (post_prec <= 0) {
        bias_cur[m] <- 0
      } else {
        post_mean <- (resid_sum / eta_m) / post_prec
        bias_cur[m] <- rnorm(1, post_mean, sqrt(1 / post_prec))
      }
    }
    bias_cur
  }
  
  # ===== MAIN MCMC LOOP =====
  # Record K_occ every iteration to make progress averaging bulletproof
  K_hist <- rep(NA_real_, n_iter)
  
  for (iter in 1:n_iter) {
    # ===== ITERATION START: WARMUP HANDLING =====
    if (iter == warmup_iters + 1 && length(revealed_idx) && warmup_iters > 0 && unpin) {
      message("Warm-up ended: released revealed curves for regular DP assignments.")
    }
    
    # ===== SLICE SAMPLING: STICK-BREAKING AND CLUSTER EXPANSION =====
    # Walker/Kalli slice expansion
    pi <- stick_to_pi(v)
    u  <- sapply(1:N, function(i) runif(1, 0, pi[z[i]]))
    u_star <- min(u)
    v <- extend_sticks_until(v, alpha, u_star)
    pi <- stick_to_pi(v)
    K  <- length(v)
    while (length(params) < K) {
      params[[length(params) + 1]] <- draw_new_cluster_params(M, P, t, kernels, wf, J, boundary)
      k_new <- length(params)
      params <- ensure_complete_cache(params, kernels, k_new, t, M)
    }
    
    # ===== CLUSTER ASSIGNMENT FUNCTION =====
    assign_one <- function(i) {
      if (length(revealed_idx) && (i %in% revealed_idx) &&
          ( !unpin || (warmup_iters > 0 && iter <= warmup_iters) )) {
        return(1L)
      }
      S <- which(pi > u[i]); if (length(S) == 0) S <- 1L
      logw <- rep(-Inf, length(S))
      for (ss in seq_along(S)) {
        k <- S[ss]
        if (is.null(params[[k]]$mu_cached) || is.null(params[[k]]$mu_cached_iter) ||
            params[[k]]$mu_cached_iter != iter) {
          mu_wave <- if (!is.null(params[[k]]$beta_ch) && length(params[[k]]$beta_ch)) {
            compute_mu_from_beta(params[[k]]$beta_ch, wf, J, boundary, P)
          } else matrix(0, P, M)
          mu_k <- add_bias_to_mu(mu_wave, params[[k]]$bias)
          # Only convert to matrix if not already a matrix or if storage mode is wrong
          if (!is.matrix(mu_k) || storage.mode(mu_k) != "double") {
            params[[k]]$mu_cached <- as_num_mat(mu_k)
          } else {
            params[[k]]$mu_cached <- mu_k
          }
          params[[k]]$mu_cached_iter <- iter
        }
        ll <- ll_curve_k(k, Y[[i]], params[[k]]$mu_cached)
        logw[ss] <- log(pi[k]) + ll
      }
      w <- exp(logw - max(logw)); w <- w / sum(w)
      S[sample.int(length(S), 1, prob = w)]
    }
    
    # ===== CLUSTER ASSIGNMENT UPDATES =====
    for (kk in seq_along(params)) params <- ensure_complete_cache(params, kernels, kk, t, M)
    if (isTRUE(diagnostics)) z_prev <- z

    for (i in 1:N) z[i] <- assign_one(i)
    
    # ===== STICK-BREAKING PARAMETER UPDATES =====
    v  <- update_v_given_z(v, z, alpha)
    pi <- stick_to_pi(v)
    K  <- length(v)
    
    # ===== CLUSTER-SPECIFIC UPDATES =====
    for (k in 1:K) {
      idx <- which(z == k)
      if (!length(idx)) next
      # ===== BIAS TERM SAMPLING (if mean_intercept=TRUE) =====
      if (isTRUE(mean_intercept)) {
        mu_wave_cur <- compute_mu_from_beta(params[[k]]$beta_ch, wf, J, boundary, P)
        params[[k]]$bias <- sample_bias(k, idx, mu_wave_cur)
      }
      # ===== WAVELET BLOCK UPDATES (BESOV PRIORS) =====
      bias_coeff <- if (isTRUE(mean_intercept)) bias_to_wavelet(params[[k]]$bias) else NULL
      upd <- update_cluster_wavelet_params_besov(
        idx = idx, precomp = precomp_all, M = M,
        wpar = params[[k]]$wpar,
        sigma2_m = params[[k]]$sigma2,
        tau_sigma = params[[k]]$tau_sigma,
        kappa_pi = kappa_pi, c2 = c2, tau_pi = tau_pi,
        g_hyp = params[[k]]$g_hyp,
        a_sig = a_sig, b_sig = b_sig, a_tau = a_tau, b_tau = b_tau,
        bias_coeff = bias_coeff
      )
      params[[k]]$wpar       <- upd$wpar
      params[[k]]$beta_ch    <- upd$beta_ch
      params[[k]]$sigma2     <- upd$sigma2_m
      params[[k]]$tau_sigma  <- upd$tau_sigma

      # ===== MEAN FUNCTION RECOMPUTATION AND CACHING =====
      mu_wave <- compute_mu_from_beta(params[[k]]$beta_ch, wf, J, boundary, P)
      if (isTRUE(mean_intercept)) {
        params[[k]]$bias <- sample_bias(k, idx, mu_wave)
      }
      mu_k <- add_bias_to_mu(mu_wave, params[[k]]$bias)
      # Only convert to matrix if not already a matrix or if storage mode is wrong
      if (!is.matrix(mu_k) || storage.mode(mu_k) != "double") {
        params[[k]]$mu_cached <- as_num_mat(mu_k)
      } else {
        params[[k]]$mu_cached <- mu_k
      }
      params[[k]]$mu_cached_iter <- iter

      # ===== KERNEL AND ICM UPDATES =====
      Yk <- lapply(idx, function(ii) Y[[ii]])
      params <- cc_switch_kernel_eig(k, params, kernels, t, Yk)  # Carlin-Chib kernel swapping
      params <- mh_update_kernel_eig(k, params, kernels, t, Yk, a_eta, b_eta)  # Kernel parameters
      params <- mh_update_L_eig(     k, params, kernels, t, Yk, mh_step_L)      # L matrix
      params <- mh_update_eta_eig(   k, params, kernels, t, Yk, mh_step_eta, a_eta, b_eta)  # eta vector
      params <- mh_update_tauB_eig(  k, params, kernels, t, Yk, mh_step_tauB)  # tau_B scalar
    }
    
    # ===== DIAGNOSTICS AND TRACKING =====
    # Occupied clusters this iteration; record to K_hist *every* iter
    Kocc <- length(unique(z))
    if (isTRUE(diagnostics)) {
      diag$ari_all[iter] <- adj_rand_index(z_prev, z)
      diag$global$K_occ_all[iter] <- Kocc
    }
    K_hist[iter] <- as.numeric(Kocc)
    
    # ===== DP CONCENTRATION PARAMETER UPDATE (ESCÓBAR-WEST) =====
    # Escobar–West update for alpha
    eta_aux <- rbeta(1, alpha + 1, N)
    mix <- (alpha_prior[1] + Kocc - 1) /
      (N * (alpha_prior[2] - log(eta_aux)) + alpha_prior[1] + Kocc - 1)
    if (runif(1) < mix) {
      alpha <- rgamma(1, alpha_prior[1] + Kocc, alpha_prior[2] - log(eta_aux))
    } else {
      alpha <- rgamma(1, alpha_prior[1] + Kocc - 1, alpha_prior[2] - log(eta_aux))
    }
    
    # ===== SAMPLE STORAGE (POST-BURN) =====
    # Save draw (if any)
    if (keep > 0 && iter > burn && ((iter - burn) %% thin == 0)) {
      sidx <- sidx + 1L
      Z_s[sidx, ] <- z
      alpha_s[sidx] <- alpha
      tab <- sort(table(z), decreasing = TRUE)
      k_big <- as.integer(names(tab)[1])
      k_sec <- if (length(tab) >= 2) as.integer(names(tab)[2]) else k_big
      kern_s[sidx] <- params[[k_big]]$kern_idx
      K_s[sidx] <- Kocc
      
      totll <- 0
      for (i in 1:N) {
        ki <- z[i]
        totll <- totll + ll_curve_k(ki, Y[[i]], params[[ki]]$mu_cached)
      }
      loglik_s[sidx] <- totll
      
      if (isTRUE(diagnostics)) {
        diag$global$K_occ[sidx]  <- Kocc
        diag$global$alpha[sidx]  <- alpha
        diag$global$loglik[sidx] <- totll
        
        snap_k <- function(k) {
          kc <- params[[k]]
          kp <- kc$thetas[[ kc$kern_idx ]]; kern_pars <- unlist(kp); names(kern_pars) <- names(kp)
          Lsel <- c(L11 = kc$L[1,1],
                    L21 = ifelse(nrow(kc$L)>=2, kc$L[2,1], NA_real_),
                    L31 = ifelse(nrow(kc$L)>=3, kc$L[3,1], NA_real_))
          pi_g <- lapply(diag$meta$monitor_levels, function(lev)
            c(pi = kc$wpar$pi_level[[lev]] %||% NA_real_,
              g  = kc$wpar$g_level[[lev]] %||% NA_real_))
          pi_g <- do.call(rbind, pi_g); rownames(pi_g) <- diag$meta$monitor_levels
          
          spars <- matrix(NA_real_, nrow=length(kc$wpar$lev_names), ncol=length(kc$wpar$gamma_ch))
          rownames(spars) <- kc$wpar$lev_names
          colnames(spars) <- paste0("ch", seq_along(kc$wpar$gamma_ch))
          zero_fw <- wt_forward_mat(matrix(0, P, M), wf, J, boundary)
          for (m in seq_along(kc$wpar$gamma_ch)) {
            gvec <- kc$wpar$gamma_ch[[m]]
            off <- 0
            for (lev in kc$wpar$lev_names) {
              ids <- zero_fw[[m]]$map$idx[[lev]]
              Llev <- length(ids)
              if (!is.null(Llev) && Llev>0) {
                idx <- (off+1):(off+Llev); off <- off+Llev
                spars[lev, m] <- mean(gvec[idx]==1)
              }
            }
          }
          list(kern_idx=kc$kern_idx, kern_pars=kern_pars, tau_B=kc$tau_B,
               Lsel=Lsel, eta=kc$eta, sigma2=kc$sigma2, tau_sigma=kc$tau_sigma,
               pi_g=pi_g, sparsity=spars, bias=kc$bias, acc=kc$acc)
        }
        diag$clusters$A[[sidx]] <- snap_k(k_big)
        diag$clusters$B[[sidx]] <- snap_k(k_sec)
      }
    }
    
    # ===== PROGRESS REPORTING =====
    # Progress line (robust; never "-" after burn)
    if (iter %% 10 == 0) {
      avg_saved <- if (keep > 0 && sidx > 0) mean(K_s[seq_len(sidx)], na.rm = TRUE) else NA_real_
      avg_iter  <- if (iter > burn) mean(K_hist[(burn + 1):iter], na.rm = TRUE) else NA_real_
      show_avg  <- if (!is.na(avg_saved)) {
        sprintf("%.2f", avg_saved)
      } else if (!is.na(avg_iter)) {
        sprintf("%.2f", avg_iter)
      } else {
        "-"
      }
      cat(sprintf(
        "\rIter %d/%d | K_occ=%d | K_sticks=%d | avgK_post=%s   ",
        iter, n_iter, Kocc, length(v), show_avg
      ))
      flush.console()
    }
  }
  # ===== MCMC LOOP COMPLETE =====
  cat("\nDone.\n")
  
  # ===== POST-PROCESSING: ARI COMPUTATION =====
  # ARI across consecutive saved draws (if any)
  if (isTRUE(diagnostics) && keep >= 2) {
    for (s in 2:keep) {
      diag$ari[s-1] <- adj_rand_index(Z_s[s-1,], Z_s[s,])
    }
  }
  
  # ===== FINAL SUMMARIES AND OUTPUT =====
  # Final summaries
  if (keep > 0 && length(K_s) > 0 && any(is.finite(K_s))) {
    cat(sprintf("Post-burn avg K (saved draws): %.3f\n", mean(K_s[is.finite(K_s)])))
  } else {
    cat("No saved draws to average over (check burn/thin/n_iter).\n")
  }
  if (any(!is.na(K_hist[(burn+1):n_iter]))) {
    cat(sprintf("Post-burn avg K (per-iter): %.3f\n",
                mean(K_hist[(burn+1):n_iter], na.rm = TRUE)))
  }
  
  # ===== RETURN RESULTS =====
  list(
    Z = Z_s, alpha = alpha_s, kern = kern_s,
    params = params, v = v, pi = stick_to_pi(v),
    revealed_idx = revealed_idx,
    K_occ = K_s, loglik = loglik_s,
    diagnostics = diag
  )
}
