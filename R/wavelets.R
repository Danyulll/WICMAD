
# Wavelet forward/inverse & precomputation -------------------------------

wt_forward_1d <- function(y, wf = "la8", J = NULL, boundary = "periodic") {
  P <- length(y); J <- ensure_dyadic_J(P, J)
  w <- waveslim::dwt(y, wf=wf, n.levels=J, boundary=boundary)
  vec <- c(w$d1); idx <- list(d1 = seq_along(w$d1)); off <- length(w$d1)
  if (J >= 2) for (lev in 2:J) {
    nm <- paste0("d", lev); v <- w[[nm]]; vec <- c(vec, v)
    idx[[nm]] <- (off + 1):(off + length(v)); off <- off + length(v)
  }
  s_nm <- paste0("s", J); vec <- c(vec, w[[s_nm]])
  idx[[s_nm]] <- (off + 1):(off + length(w[[s_nm]]))
  list(coeff = as.numeric(vec), map = list(J=J, wf=wf, boundary=boundary, P=P, idx=idx))
}
wt_inverse_1d <- function(coeff_vec, map) {
  J <- map$J
  w <- vector("list", J + 1L)
  names(w) <- c(paste0("d", 1:J), paste0("s", J))
  for (lev in 1:J) {
    nm <- paste0("d", lev); ids <- map$idx[[nm]]
    w[[nm]] <- as.numeric(coeff_vec[ids])
  }
  ids_s <- map$idx[[paste0("s", J)]]; w[[paste0("s", J)]] <- as.numeric(coeff_vec[ids_s])
  attr(w, "wavelet") <- map$wf; attr(w, "boundary") <- map$boundary; class(w) <- "dwt"
  waveslim::idwt(w)
}
wt_forward_mat <- function(y_mat, wf="la8", J=NULL, boundary="periodic") {
  M <- ncol(y_mat); out <- vector("list", M)
  for (m in 1:M) out[[m]] <- wt_forward_1d(y_mat[,m], wf, J, boundary)
  out
}
precompute_wavelets <- function(Y_list, wf, J, boundary) {
  lapply(Y_list, function(mat) wt_forward_mat(mat, wf, J, boundary))
}
stack_D_from_precomp <- function(precomp, idx, M) {
  ncoeff <- length(precomp[[ idx[1] ]][[1]]$coeff)
  N <- length(idx)
  D_arr <- array(NA_real_, dim=c(ncoeff, N, M))
  for (jj in seq_along(idx)) {
    i <- idx[jj]
    for (m in 1:M) D_arr[, jj, m] <- precomp[[i]][[m]]$coeff
  }
  maps <- precomp[[ idx[1] ]]
  list(D_arr = D_arr, maps = maps)
}
compute_mu_from_beta <- function(beta_ch, wf, J, boundary, P) {
  M <- length(beta_ch); zeros <- matrix(0, nrow=P, ncol=M)
  tmpl <- wt_forward_mat(zeros, wf, J, boundary); mu <- matrix(0, P, M)
  for (m in 1:M) mu[,m] <- wt_inverse_1d(beta_ch[[m]], tmpl[[m]]$map)
  mu
}
