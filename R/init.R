
# Initialization helpers -------------------------------------------------

.draw_empty_acc <- function(M, kernels) {
  out <- list()
  out$kernel <- list()
  pnames <- unique(unlist(lapply(kernels, function(kc) kc$pnames)))
  for (pn in pnames) out$kernel[[pn]] <- c(a=0, n=0)
  out$L    <- c(a=0, n=0)
  out$eta  <- cbind(a=rep(0,M), n=rep(0,M))
  out$tauB <- c(a=0, n=0)
  out
}

draw_new_cluster_params <- function(M, P, t, kernels, wf="la8", J=NULL, boundary="periodic") {
  if (is.null(J)) J <- floor(log2(P))
  zeros <- matrix(0, nrow=P, ncol=M)
  tmp   <- wt_forward_mat(zeros, wf=wf, J=J, boundary=boundary)
  lev_names <- names(tmp[[1]]$map$idx)
  det_names <- lev_names[grepl("^d", lev_names)]
  ncoeff    <- length(tmp[[1]]$coeff)
  beta_ch <- lapply(1:M, function(m) numeric(ncoeff))
  pi_level <- setNames(rep(0.5, length(det_names)), det_names)
  g_level  <- setNames(rep(2.0, length(det_names)), det_names)
  gamma_ch <- lapply(1:M, function(m) rbinom(ncoeff, 1, 0.2))
  thetas   <- lapply(kernels, function(kc) kc$pstar())
  list(
    wpar = list(lev_names=lev_names, pi_level=pi_level, g_level=g_level, gamma_ch=gamma_ch),
    kern_idx = sample.int(length(kernels),1),
    thetas   = thetas,
    L = diag(M),
    eta = rep(0.05, M),
    tau_B = 1.0,
    beta_ch = beta_ch,
    sigma2  = rep(0.05, M),
    tau_sigma = 1.0,
    cache = list(U=NULL, lam=NULL, V=NULL, s=NULL, last_kx_hash=""),
    mu_cached = NULL, mu_cached_iter = NA_integer_,
    acc = .draw_empty_acc(M, kernels)
  )
}

ensure_complete_cache <- function(params, kernels, k, t, M) {
  stopifnot(k >= 1L, k <= length(params))
  if (!is.numeric(params[[k]]$eta) || length(params[[k]]$eta) != M ||
      any(!is.finite(params[[k]]$eta))) {
    params[[k]]$eta <- rep(0.05, M)
  }
  if (!is.numeric(params[[k]]$tau_B) || length(params[[k]]$tau_B) != 1L ||
      !is.finite(params[[k]]$tau_B) || params[[k]]$tau_B <= 0) {
    params[[k]]$tau_B <- 1.0
  }
  if (is.null(params[[k]]$kern_idx) ||
      params[[k]]$kern_idx < 1L || params[[k]]$kern_idx > length(kernels)) {
    params[[k]]$kern_idx <- 1L
  }
  if (is.null(params[[k]]$thetas) || length(params[[k]]$thetas) != length(kernels)) {
    params[[k]]$thetas <- lapply(kernels, function(kc) kc$pstar())
  }
  if (is.null(params[[k]]$thetas[[ params[[k]]$kern_idx ]])) {
    params[[k]]$thetas[[ params[[k]]$kern_idx ]] <- kernels[[ params[[k]]$kern_idx ]]$pstar()
  }

  eeB <- eig_Bshape(params[[k]]$L, M)
  params[[k]]$cache$B_U   <- as.matrix(eeB$U)
  params[[k]]$cache$B_lam <- as.numeric(eeB$lam)

  kc <- kernels[[ params[[k]]$kern_idx ]]
  kp <- params[[k]]$thetas[[ params[[k]]$kern_idx ]]
  Kx <- kc$fun(t, kp)
  ekx <- eig_Kx(Kx)
  params[[k]]$cache$Ux_x <- as.matrix(ekx$V)
  params[[k]]$cache$lam_x <- as.numeric(ekx$s)
  vals <- unlist(kp); keys <- names(kp); if (is.null(keys)) keys <- sprintf("p%d", seq_along(vals))
  params[[k]]$cache$last_kx_hash <- paste(sprintf("%s=%.10f", keys, vals), collapse = "|")
  params
}
