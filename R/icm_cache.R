
# Fast ICM likelihood via block cache -----------------------------------

.build_icm_cache <- function(t, kern_cfg, kp, L, eta, tau_B, cache = NULL) {
  key_kx <- paste(kern_cfg$name, paste(unlist(kp), collapse="|"), sep="::")
  key_B  <- paste(round(L,8), collapse="|")
  need_kx <- is.null(cache$key_kx) || cache$key_kx != key_kx
  need_B  <- is.null(cache$key_B)  || cache$key_B  != key_B || is.null(cache$Bshape)
  need_tau_eta <- is.null(cache$tau) || cache$tau != tau_B || is.null(cache$eta) || any(cache$eta != eta)

  if (need_kx) {
    Kx <- kern_cfg$fun(t, kp)
    P_exp <- nloc(t)
    if (!is.matrix(Kx) || nrow(Kx) != P_exp || ncol(Kx) != P_exp) {
      stop(sprintf("Kernel Kx wrong size: got %dx%d, expected %dx%d.", nrow(Kx), ncol(Kx), P_exp, P_exp))
    }
    eig <- eigen(Kx, symmetric = TRUE)
    cache$Ux_x   <- eig$vectors
    cache$lam_x  <- pmax(eig$values, 1e-12)
    cache$key_kx <- key_kx
  }

  if (need_B) {
    Bshape <- tcrossprod(L)
    trB <- sum(diag(Bshape))
    if (trB > 0) Bshape <- Bshape * (nrow(Bshape) / trB)
    cache$Bshape <- Bshape
    cache$key_B  <- key_B
  }

  if (need_kx || need_B || need_tau_eta || is.null(cache$chol_list)) {
    if (is.null(cache$lam_x)) stop("Cache missing lam_x.")
    P_need <- length(cache$lam_x)
    M <- length(eta)
    cache$chol_list <- vector("list", P_need)
    logdet_sum <- 0
    Deta <- diag(eta, M, M)
    for (j in 1:P_need) {
      Sj <- tau_B * cache$lam_x[j] * cache$Bshape + Deta
      ok <- FALSE; tries <- 0; jitter <- 1e-8; Lj <- NULL
      while(!ok && tries < 6) {
        tries <- tries + 1
        out <- try(chol(Sj + diag(jitter, M)), silent=TRUE)
        if (!inherits(out, "try-error")) { Lj <- out; ok <- TRUE } else jitter <- jitter*10
      }
      if (!ok) stop("Cholesky failed in Σ_j block.")
      cache$chol_list[[j]] <- Lj
      logdet_sum <- logdet_sum + 2*sum(log(diag(Lj)))
    }
    cache$logdet_sum <- logdet_sum
    cache$tau <- tau_B
    cache$eta <- eta
  }
  cache
}

fast_icm_loglik_curve <- function(y_resid, cache) {
  if (is.null(cache$Ux_x) || is.null(cache$chol_list))
    stop("Cache not built: missing Ux_x or chol_list.")
  P_y  <- nrow(y_resid)
  P_ch <- length(cache$chol_list)
  if (P_y != P_ch)
    stop(sprintf("Cache dimension mismatch: y rows=%d, chol blocks=%d.", P_y, P_ch))

  Ytil <- t(cache$Ux_x) %*% y_resid
  quad <- 0
  for (j in 1:P_ch) {
    Lj <- cache$chol_list[[j]]
    v  <- as.numeric(Ytil[j, ])
    w  <- backsolve(Lj, forwardsolve(t(Lj), v, upper.tri=TRUE, transpose=TRUE))
    quad <- quad + sum(w*w)
  }
  -0.5 * (P_ch * ncol(y_resid) * log(2*base::pi) + cache$logdet_sum + quad)
}
