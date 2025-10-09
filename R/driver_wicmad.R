#' WICMAD main model (DP–ICM–GP with Besov wavelet block)
#' @param mean_intercept logical; if `TRUE` (default) include and sample a per-channel intercept
#'   in each cluster while disabling bias-augmented kernel variants.
#' @export
wicmad <- function(
    Y, t,
    n_iter = 6000, burn = 3000, thin = 5,
    alpha_prior = c(10, 1),
    wf = "la8", J = NULL, boundary = "periodic",
    mh_step_L = 0.03, mh_step_eta = 0.10, mh_step_tauB = 0.15,
    revealed_idx = integer(0),
    K_init = 5,
    use_parallel = FALSE,
    warmup_iters = 100,
    unpin = TRUE,
    kappa_pi = 0.6, c2 = 1.0, tau_pi = 40,
    a_sig = 2.5, b_sig = 0.02, a_tau = 2.0, b_tau = 2.0,
    a_eta = 2, b_eta = 0.1,
    diagnostics = TRUE,
    track_ids = NULL,
    monitor_levels = NULL,
    mean_intercept = TRUE
) {
  suppressPackageStartupMessages({
    requireNamespace("waveslim")
    requireNamespace("mvtnorm")
    requireNamespace("MASS")
    requireNamespace("invgamma")
  })
  
  # ---- Guards so we always have post-burn iterations and (if feasible) a saved draw
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
  
  N <- length(Y); P <- nrow(Y[[1]]); M <- ncol(Y[[1]]); J <- ensure_dyadic_J(P, J)
  t <- normalize_t(t, P); t <- scale_t01(t)
  
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
  
  kernels <- make_kernels(add_bias_variants = !isTRUE(mean_intercept))
  alpha <- rgamma(1, shape = alpha_prior[1], rate = alpha_prior[2])
  v     <- rbeta(K_init, 1, alpha)
  pi    <- stick_to_pi(v)
  K     <- length(v)
  
  params <- vector("list", K)
  for (k in 1:K) {
    params[[k]] <- draw_new_cluster_params(M, P, t, kernels, wf, J, boundary)
    params <- ensure_complete_cache(params, kernels, k, t, M)
  }
  
  z <- sample.int(K, N, replace = TRUE)
  if (length(revealed_idx)) z[revealed_idx] <- 1L
  
  keep <- max(0L, floor((n_iter - burn) / thin))
  Z_s <- if (keep > 0) matrix(NA_integer_, keep, N) else matrix(NA_integer_, 0, N)
  alpha_s <- if (keep > 0) numeric(keep) else numeric(0)
  kern_s  <- if (keep > 0) integer(keep) else integer(0)
  K_s     <- if (keep > 0) integer(keep) else integer(0)
  loglik_s<- if (keep > 0) numeric(keep) else numeric(0)
  
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
                    loglik=if (keep>0) rep(NA_real_, keep) else numeric(0)),
      clusters = list(A=list(), B=list()),
      ari = if (keep>1) rep(NA_real_, keep-1) else numeric(0)
    )
  }
  
  sidx <- 0L
  Y <- lapply(Y, as_num_mat)
  precomp_all <- precompute_wavelets(Y, wf, J, boundary)
  message(sprintf("Starting MCMC (N=%d, P=%d, M=%d, iters=%d; burn=%d, thin=%d, keep=%d)",
                  N, P, M, n_iter, burn, thin, keep))
  
  add_bias_to_mu <- function(mu_wave, bias_vec) {
    if (!isTRUE(mean_intercept)) return(mu_wave)
    if (length(bias_vec) == 0) return(mu_wave)
    mu_wave + matrix(rep(bias_vec, each = nrow(mu_wave)), nrow(mu_wave), length(bias_vec))
  }

  bias_to_wavelet <- function(bias_vec) {
    if (!isTRUE(mean_intercept)) return(NULL)
    if (length(bias_vec) == 0L) return(NULL)
    if (!any(abs(bias_vec) > 1e-12)) return(NULL)
    bias_mat <- matrix(rep(bias_vec, each = P), nrow = P, ncol = M)
    fw <- wt_forward_mat(bias_mat, wf, J, boundary)
    lapply(fw, function(comp) comp$coeff)
  }

  ll_curve_k <- function(k, Yi, mu_k) {
    kc <- kernels[[ params[[k]]$kern_idx ]]
    kp <- params[[k]]$thetas[[ params[[k]]$kern_idx ]]
    cache <- .build_icm_cache(t, kc, kp, params[[k]]$L, params[[k]]$eta, params[[k]]$tau_B, params[[k]]$cache)
    params[[k]]$cache <<- cache
    y_resid <- Yi - mu_k
    fast_icm_loglik_curve(y_resid, cache)
  }

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
  
  # ---- Record K_occ every iteration to make progress averaging bulletproof
  K_hist <- rep(NA_real_, n_iter)
  
  for (iter in 1:n_iter) {
    if (iter == warmup_iters + 1 && length(revealed_idx) && warmup_iters > 0 && unpin) {
      message("Warm-up ended: released revealed curves for regular DP assignments.")
    }
    
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
          params[[k]]$mu_cached <- as_num_mat(mu_k)
          params[[k]]$mu_cached_iter <- iter
        }
        ll <- ll_curve_k(k, Y[[i]], params[[k]]$mu_cached)
        logw[ss] <- log(pi[k]) + ll
      }
      w <- exp(logw - max(logw)); w <- w / sum(w)
      S[sample.int(length(S), 1, prob = w)]
    }
    
    for (kk in seq_along(params)) params <- ensure_complete_cache(params, kernels, kk, t, M)
    for (i in 1:N) z[i] <- assign_one(i)
    
    v  <- update_v_given_z(v, z, alpha)
    pi <- stick_to_pi(v)
    K  <- length(v)
    
    for (k in 1:K) {
      idx <- which(z == k)
      if (!length(idx)) next
      if (isTRUE(mean_intercept)) {
        mu_wave_cur <- compute_mu_from_beta(params[[k]]$beta_ch, wf, J, boundary, P)
        params[[k]]$bias <- sample_bias(k, idx, mu_wave_cur)
      }
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

      mu_wave <- compute_mu_from_beta(params[[k]]$beta_ch, wf, J, boundary, P)
      if (isTRUE(mean_intercept)) {
        params[[k]]$bias <- sample_bias(k, idx, mu_wave)
      }
      mu_k <- add_bias_to_mu(mu_wave, params[[k]]$bias)
      params[[k]]$mu_cached      <- as_num_mat(mu_k)
      params[[k]]$mu_cached_iter <- iter

      Yk <- lapply(idx, function(ii) Y[[ii]])
      params <- cc_switch_kernel_eig(k, params, kernels, t, Yk)
      params <- mh_update_kernel_eig(k, params, kernels, t, Yk, a_eta, b_eta)
      params <- mh_update_L_eig(     k, params, kernels, t, Yk, mh_step_L)
      params <- mh_update_eta_eig(   k, params, kernels, t, Yk, mh_step_eta, a_eta, b_eta)
      params <- mh_update_tauB_eig(  k, params, kernels, t, Yk, mh_step_tauB)
    }
    
    # ---- Occupied clusters this iteration; record to K_hist *every* iter
    Kocc <- length(unique(z))
    K_hist[iter] <- as.numeric(Kocc)
    
    # ---- Escobar–West update for alpha
    eta_aux <- rbeta(1, alpha + 1, N)
    mix <- (alpha_prior[1] + Kocc - 1) /
      (N * (alpha_prior[2] - log(eta_aux)) + alpha_prior[1] + Kocc - 1)
    if (runif(1) < mix) {
      alpha <- rgamma(1, alpha_prior[1] + Kocc, alpha_prior[2] - log(eta_aux))
    } else {
      alpha <- rgamma(1, alpha_prior[1] + Kocc - 1, alpha_prior[2] - log(eta_aux))
    }
    
    # ---- Save draw (if any)
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
    
    # ---- Progress line (robust; never "-" after burn)
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
  cat("\nDone.\n")
  
  # ---- ARI across consecutive saved draws (if any)
  if (isTRUE(diagnostics) && keep >= 2) {
    adj_rand_index <- function(z1, z2) {
      tab <- table(z1, z2)
      sum_comb <- sum(choose(tab, 2))
      a <- rowSums(tab); b <- colSums(tab)
      sum_a <- sum(choose(a, 2)); sum_b <- sum(choose(b, 2))
      tot <- choose(length(z1), 2)
      (sum_comb - (sum_a*sum_b)/tot) / (0.5*(sum_a+sum_b) - (sum_a*sum_b)/tot)
    }
    for (s in 2:keep) {
      diag$ari[s-1] <- adj_rand_index(Z_s[s-1,], Z_s[s,])
    }
  }
  
  # ---- Final summaries
  if (keep > 0 && length(K_s) > 0 && any(is.finite(K_s))) {
    cat(sprintf("Post-burn avg K (saved draws): %.3f\n", mean(K_s[is.finite(K_s)])))
  } else {
    cat("No saved draws to average over (check burn/thin/n_iter).\n")
  }
  if (any(!is.na(K_hist[(burn+1):n_iter]))) {
    cat(sprintf("Post-burn avg K (per-iter): %.3f\n",
                mean(K_hist[(burn+1):n_iter], na.rm = TRUE)))
  }
  
  list(
    Z = Z_s, alpha = alpha_s, kern = kern_s,
    params = params, v = v, pi = stick_to_pi(v),
    revealed_idx = revealed_idx,
    K_occ = K_s, loglik = loglik_s,
    diagnostics = diag
  )
}
