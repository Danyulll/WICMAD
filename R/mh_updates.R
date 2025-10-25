
# MCMC update steps ------------------------------------------------------

#' Sum Log-Likelihoods Across Curves
#'
#' Computes the total log-likelihood across multiple curves using ICM cache.
#'
#' @param curves List of residual curve matrices
#' @param cache ICM cache object
#' @return Numeric total log-likelihood
sum_ll_curves <- function(curves, cache) {
  # Use C++ batch implementation for maximum performance
  if (requireNamespace("RcppEigen", quietly = TRUE)) {
    tryCatch({
      return(sum(fast_icm_loglik_curves_batch(curves, cache$Ux_x, cache$chol_list, cache$logdet_sum)))
    }, error = function(e) {
      warning("C++ batch implementation failed, falling back to R: ", e$message)
    })
  }
  
  # Fallback to optimized R implementation
  sum(sapply(curves, function(y) fast_icm_loglik_curve(y, cache)))
}

#' Update Kernel Parameters via Metropolis-Hastings
#'
#' Performs Metropolis-Hastings updates for kernel hyperparameters using eigendecomposition.
#'
#' @param k Integer cluster index
#' @param params List of cluster parameters
#' @param kernels List of kernel configurations
#' @param t Numeric vector or matrix of input coordinates
#' @param Y_list List of observed curve matrices
#' @param a_eta Numeric shape parameter for eta prior
#' @param b_eta Numeric rate parameter for eta prior
#' @return Updated params list
mh_update_kernel_eig <- function(k, params, kernels, t, Y_list, a_eta, b_eta) {
  kc <- kernels[[ params[[k]]$kern_idx ]]
  kp <- params[[k]]$thetas[[ params[[k]]$kern_idx ]]
  curves <- lapply(Y_list, function(Yi) Yi - params[[k]]$mu_cached)
  for (pn in kc$pnames) {
    cur <- kp[[pn]]
    if (pn == "period") {
      z_cur <- logit(pmin(pmax(cur, 1e-6), 1 - 1e-6))
      z_prp <- rnorm(1, z_cur, kc$prop_sd[[pn]]); prp <- invlogit(z_prp)
      kp_prop <- kp; kp_prop[[pn]] <- prp
      cache_cur <- .build_icm_cache(t, kc, kp, params[[k]]$L, params[[k]]$eta, params[[k]]$tau_B, params[[k]]$cache)
      ll_cur <- sum_ll_curves(curves, cache_cur)
      cache_prp <- .build_icm_cache(t, kc, kp_prop, params[[k]]$L, params[[k]]$eta, params[[k]]$tau_B, params[[k]]$cache)
      ll_prp <- sum_ll_curves(curves, cache_prp)
      lp_cur <- kc$prior(kp); lp_prp <- kc$prior(kp_prop)
      a <- (ll_prp + lp_prp) - (ll_cur + lp_cur)
      if (is.finite(a) && log(runif(1)) < a) {
        kp <- kp_prop; params[[k]]$thetas[[ params[[k]]$kern_idx ]] <- kp; params[[k]]$cache <- cache_prp
        params[[k]]$acc$kernel[[pn]]["a"] <- params[[k]]$acc$kernel[[pn]]["a"] + 1
      }
      params[[k]]$acc$kernel[[pn]]["n"] <- params[[k]]$acc$kernel[[pn]]["n"] + 1
    } else {
      prp <- rlnorm(1, log(cur), kc$prop_sd[[pn]]); kp_prop <- kp; kp_prop[[pn]] <- prp
      cache_cur <- .build_icm_cache(t, kc, kp, params[[k]]$L, params[[k]]$eta, params[[k]]$tau_B, params[[k]]$cache)
      ll_cur <- sum_ll_curves(curves, cache_cur)
      cache_prp <- .build_icm_cache(t, kc, kp_prop, params[[k]]$L, params[[k]]$eta, params[[k]]$tau_B, params[[k]]$cache)
      ll_prp <- sum_ll_curves(curves, cache_prp)
      lp_cur <- kc$prior(kp); lp_prp <- kc$prior(kp_prop)
      q_cgpr <- dlnorm(cur, log(prp), kc$prop_sd[[pn]], log=TRUE)
      q_prgc <- dlnorm(prp, log(cur), kc$prop_sd[[pn]], log=TRUE)
      a <- (ll_prp + lp_prp + q_cgpr) - (ll_cur + lp_cur + q_prgc)
      if (is.finite(a) && log(runif(1)) < a) {
        kp <- kp_prop; params[[k]]$thetas[[ params[[k]]$kern_idx ]] <- kp; params[[k]]$cache <- cache_prp
        params[[k]]$acc$kernel[[pn]]["a"] <- params[[k]]$acc$kernel[[pn]]["a"] + 1
      }
      params[[k]]$acc$kernel[[pn]]["n"] <- params[[k]]$acc$kernel[[pn]]["n"] + 1
    }
  }
  params
}

mh_update_L_eig <- function(k, params, kernels, t, Y_list, mh_step_L) {
  th  <- pack_L(params[[k]]$L)
  thp <- th + rnorm(length(th), 0, mh_step_L)
  Lp  <- unpack_L(thp, nrow(params[[k]]$L))
  curves <- lapply(Y_list, function(Yi) Yi - params[[k]]$mu_cached)
  kc <- kernels[[ params[[k]]$kern_idx ]]; kp <- params[[k]]$thetas[[ params[[k]]$kern_idx ]]
  cache_cur <- .build_icm_cache(t, kc, kp, params[[k]]$L, params[[k]]$eta, params[[k]]$tau_B, params[[k]]$cache)
  ll_cur <- sum_ll_curves(curves, cache_cur)
  cache_prp <- .build_icm_cache(t, kc, kp, Lp,                  params[[k]]$eta, params[[k]]$tau_B, params[[k]]$cache)
  ll_prp <- sum_ll_curves(curves, cache_prp)
  lp_cur <- sum(dnorm(th,  0, 1, log=TRUE)); lp_prp <- sum(dnorm(thp, 0, 1, log=TRUE))
  a <- (ll_prp + lp_prp) - (ll_cur + lp_cur)
  if (is.finite(a) && log(runif(1)) < a) {
    params[[k]]$L <- Lp; params[[k]]$cache <- cache_prp
    params[[k]]$acc$L["a"] <- params[[k]]$acc$L["a"] + 1
  }
  params[[k]]$acc$L["n"] <- params[[k]]$acc$L["n"] + 1
  params
}

mh_update_eta_eig <- function(k, params, kernels, t, Y_list, mh_step_eta, a_eta, b_eta) {
  curves <- lapply(Y_list, function(Yi) Yi - params[[k]]$mu_cached)
  kc <- kernels[[ params[[k]]$kern_idx ]]; kp <- params[[k]]$thetas[[ params[[k]]$kern_idx ]]
  for (j in seq_along(params[[k]]$eta)) {
    cur <- params[[k]]$eta[j]
    prp <- rlnorm(1, log(cur), mh_step_eta); etap <- params[[k]]$eta; etap[j] <- prp
    cache_cur <- .build_icm_cache(t, kc, kp, params[[k]]$L, params[[k]]$eta,  params[[k]]$tau_B, params[[k]]$cache)
    cache_prp <- .build_icm_cache(t, kc, kp, params[[k]]$L, etap,               params[[k]]$tau_B, params[[k]]$cache)
    ll_cur <- sum_ll_curves(curves, cache_cur)
    ll_prp <- sum_ll_curves(curves, cache_prp)
    lp_cur <- invgamma::dinvgamma(cur, a_eta, b_eta, log=TRUE)
    lp_prp <- invgamma::dinvgamma(prp, a_eta, b_eta, log=TRUE)
    q_cgpr <- dlnorm(cur, log(prp), mh_step_eta, log=TRUE); q_prgc <- dlnorm(prp, log(cur), mh_step_eta, log=TRUE)
    a <- (ll_prp + lp_prp + q_cgpr) - (ll_cur + lp_cur + q_prgc)
    if (is.finite(a) && log(runif(1)) < a) { params[[k]]$eta <- etap; params[[k]]$cache <- cache_prp
    params[[k]]$acc$eta[j,"a"] <- params[[k]]$acc$eta[j,"a"] + 1 }
    params[[k]]$acc$eta[j,"n"] <- params[[k]]$acc$eta[j,"n"] + 1
  }
  params
}

mh_update_tauB_eig <- function(k, params, kernels, t, Y_list, mh_step_tauB) {
  curves <- lapply(Y_list, function(Yi) Yi - params[[k]]$mu_cached)
  kc <- kernels[[ params[[k]]$kern_idx ]]; kp <- params[[k]]$thetas[[ params[[k]]$kern_idx ]]
  cur <- params[[k]]$tau_B
  prp <- rlnorm(1, log(cur), mh_step_tauB)
  cache_cur <- .build_icm_cache(t, kc, kp, params[[k]]$L, params[[k]]$eta, cur, params[[k]]$cache)
  cache_prp <- .build_icm_cache(t, kc, kp, params[[k]]$L, params[[k]]$eta, prp, params[[k]]$cache)
  ll_cur <- sum_ll_curves(curves, cache_cur)
  ll_prp <- sum_ll_curves(curves, cache_prp)
  lp_cur <- -log1p(cur); lp_prp <- -log1p(prp)
  q_cgpr <- dlnorm(cur, log(prp), mh_step_tauB, log=TRUE); q_prgc <- dlnorm(prp, log(cur), mh_step_tauB, log=TRUE)
  a <- (ll_prp + lp_prp + q_cgpr) - (ll_cur + lp_cur + q_prgc)
  if (is.finite(a) && log(runif(1)) < a) { params[[k]]$tau_B <- prp; params[[k]]$cache <- cache_prp
  params[[k]]$acc$tauB["a"] <- params[[k]]$acc$tauB["a"] + 1 }
  params[[k]]$acc$tauB["n"] <- params[[k]]$acc$tauB["n"] + 1
  params
}

cc_switch_kernel_eig <- function(k, params, kernels, t, Y_list) {
  Mmod <- length(kernels)
  p_m  <- rep(1 / Mmod, Mmod)
  theta_draws <- vector("list", Mmod)
  for (m in 1:Mmod) {
    theta_draws[[m]] <- if (m == params[[k]]$kern_idx) params[[k]]$thetas[[m]] else kernels[[m]]$pstar()
  }
  curves <- lapply(Y_list, function(Yi) Yi - params[[k]]$mu_cached)
  logw <- rep(NA_real_, Mmod)
  for (m in 1:Mmod) {
    kc_m <- kernels[[m]]; kp_m <- theta_draws[[m]]
    cache_m <- .build_icm_cache(t, kc_m, kp_m, params[[k]]$L, params[[k]]$eta, params[[k]]$tau_B, params[[k]]$cache)
    ll_m <- sum_ll_curves(curves, cache_m)
    logw[m] <- log(p_m[m]) + ll_m
  }
  w <- exp(logw - max(logw)); w <- w / sum(w)
  new_idx <- sample.int(Mmod, 1, prob = w)
  params[[k]]$kern_idx <- new_idx
  params[[k]]$thetas   <- theta_draws
  kc <- kernels[[new_idx]]; kp <- params[[k]]$thetas[[new_idx]]
  params[[k]]$cache <- .build_icm_cache(t, kc, kp, params[[k]]$L, params[[k]]$eta, params[[k]]$tau_B, params[[k]]$cache)
  params
}
