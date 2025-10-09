
# Wavelet block with Besov priors ---------------------------------------

update_cluster_wavelet_params_besov <- function(
    idx, precomp, M, wpar, sigma2_m, tau_sigma,
    kappa_pi = 0.6, c2 = 1.0, tau_pi = 40,
    g_hyp = NULL,
    a_sig = 2.5, b_sig = 0.02,   # for sigma2_m
    a_tau = 2.0, b_tau = 2.0     # for tau_sigma
) {
  if (length(idx) == 0) {
    return(list(
      wpar = wpar,
      beta_ch = wpar$beta_ch %||% lapply(1:M, function(m) 0),
      sigma2_m = sigma2_m,
      tau_sigma = tau_sigma,
      maps = NULL
    ))
  }
  if (length(sigma2_m) != M) stop("sigma2_m must be length M.")
  stk <- stack_D_from_precomp(precomp, idx, M)
  D   <- stk$D_arr; maps <- stk$maps
  ncoeff <- dim(D)[1]; N <- dim(D)[2]
  lev_names <- names(maps[[1]]$map$idx)
  det_names <- lev_names[grepl("^d", lev_names)]
  s_name    <- lev_names[grepl("^s", lev_names)]
  if (is.null(wpar$pi_level))  wpar$pi_level <- setNames(rep(0.5, length(det_names)), det_names)
  if (is.null(wpar$g_level))   wpar$g_level  <- setNames(rep(2.0, length(det_names)), det_names)
  if (is.null(wpar$gamma_ch))  wpar$gamma_ch <- lapply(1:M, function(m) rbinom(ncoeff, 1, 0.2))

  # 1) gamma updates (per level & channel)
  for (m in 1:M) {
    Dm <- matrix(D[ , , m, drop=FALSE], nrow=ncoeff, ncol=N)
    gam <- wpar$gamma_ch[[m]]
    for (lev in det_names) {
      ids <- maps[[m]]$map$idx[[lev]]; if (!length(ids)) next
      pi_j <- wpar$pi_level[[lev]]; g_j  <- wpar$g_level[[lev]]
      Dsub <- Dm[ids, , drop=FALSE]

      v_spike <- sigma2_m[m]; v_slab <- (1 + g_j) * sigma2_m[m]
      ll_spike <- -0.5 * rowSums( log(2*base::pi*v_spike) + (Dsub^2)/v_spike )
      ll_slab  <- -0.5 * rowSums( log(2*base::pi*v_slab ) + (Dsub^2)/v_slab  )
      logit_val <- log(pi_j) + ll_slab - (log(1-pi_j) + ll_spike)
      p1 <- plogis(pmax(pmin(logit_val,35), -35))

      gam[ids] <- rbinom(length(ids), 1, p1)
    }
    if (length(s_name) == 1) {
      ids_s <- maps[[m]]$map$idx[[s_name]]
      if (length(ids_s) > 0) gam[ids_s] <- 1L
    }
    wpar$gamma_ch[[m]] <- gam
  }

  # 2) g_level inverse-gamma updates
  for (lev in det_names) {
    shape0 <- if (!is.null(g_hyp)) g_hyp[lev, "shape"] else 2.0
    rate0  <- if (!is.null(g_hyp)) g_hyp[lev, "rate"]  else 2.0
    ss_over_sigma <- 0; n_sel_total <- 0
    for (m in 1:M) {
      ids <- maps[[m]]$map$idx[[lev]]; if (!length(ids)) next
      Dm <- matrix(D[ , , m, drop=FALSE], nrow=ncoeff, ncol=N)
      sel <- wpar$gamma_ch[[m]][ids] == 1
      if (any(sel)) {
        Did <- Dm[ids[sel], , drop=FALSE]
        ss_over_sigma <- ss_over_sigma + sum(Did^2) / sigma2_m[m]
        n_sel_total   <- n_sel_total + nrow(Did)
      }
    }
    shape_post <- shape0 + 0.5 * n_sel_total
    rate_post  <- rate0  + 0.5 * ss_over_sigma
    wpar$g_level[[lev]] <- invgamma::rinvgamma(1, shape=shape_post, rate=rate_post)
  }

  # 3) pi_level (Besov prior)
  for (lev in det_names) {
    jnum <- as.integer(sub("^d", "", lev))
    m_j  <- max(min(kappa_pi * 2^(-c2 * jnum), 1 - 1e-6), 1e-6)
    a0   <- tau_pi * m_j
    b0   <- tau_pi * (1 - m_j)
    n1 <- 0; n0 <- 0
    for (m in 1:M) {
      ids <- maps[[m]]$map$idx[[lev]]; if (!length(ids)) next
      gm <- wpar$gamma_ch[[m]][ids]; n1 <- n1 + sum(gm == 1); n0 <- n0 + sum(gm == 0)
    }
    wpar$pi_level[[lev]] <- rbeta(1, a0 + n1, b0 + n0)
  }

  # 4) beta sampling
  beta_ch <- lapply(1:M, function(m) numeric(ncoeff))
  for (m in 1:M) {
    Dm <- matrix(D[ , , m, drop=FALSE], nrow=ncoeff, ncol=N)
    gam <- wpar$gamma_ch[[m]]; b <- numeric(ncoeff)
    for (lev in det_names) {
      ids <- maps[[m]]$map$idx[[lev]]; if (!length(ids)) next
      g_j <- wpar$g_level[[lev]]
      n   <- N
      Dbar<- rowMeans(Dm[ids, , drop=FALSE])
      shrink <- n / (n + 1 / g_j)
      mean_post <- shrink * Dbar
      var_post  <- sigma2_m[m] / (n + 1 / g_j)
      is_on <- (gam[ids] == 1L)
      if (any(is_on)) {
        k <- sum(is_on)
        b[ids[is_on]] <- rnorm(k, mean_post[is_on], sqrt(var_post))
      }
    }
    if (length(s_name) == 1) {
      ids_s <- maps[[m]]$map$idx[[s_name]]
      if (length(ids_s) > 0) {
        Dbar_s <- rowMeans(Dm[ids_s, , drop=FALSE])
        b[ids_s] <- rnorm(length(ids_s), Dbar_s, sqrt(sigma2_m[m] / N))
      }
    }
    beta_ch[[m]] <- b
  }

  # 5) sigma2_m updates + 6) tau_sigma
  sigma2_m_new <- sigma2_m
  n_eff_m <- ncoeff * N
  for (m in 1:M) {
    Dm <- matrix(D[ , , m, drop=FALSE], nrow=ncoeff, ncol=N)
    resid <- sweep(Dm, 1, beta_ch[[m]], FUN = "-")
    ss_m <- sum(resid^2)
    shape_post <- a_sig + 0.5 * n_eff_m
    rate_post  <- b_sig * tau_sigma + 0.5 * ss_m
    sigma2_m_new[m] <- invgamma::rinvgamma(1, shape = shape_post, rate = rate_post)
  }
  a_post <- a_tau + M * a_sig
  b_post <- b_tau + b_sig * sum(1 / sigma2_m_new)
  tau_sigma_new <- rgamma(1, shape = a_post, rate = b_post)

  list(wpar=wpar, beta_ch=beta_ch, sigma2_m=sigma2_m_new, tau_sigma=tau_sigma_new, maps=maps)
}
