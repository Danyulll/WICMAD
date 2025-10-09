
# Utilities ---------------------------------------------------------------

scale_t01 <- function(t) {
  if (is.vector(t)) {
    rng <- max(t) - min(t); if (rng <= 0) return(as.numeric(t))
    return((as.numeric(t) - min(t)) / rng)
  }
  tt <- as.matrix(t)
  for (j in seq_len(ncol(tt))) {
    mn <- min(tt[, j]); mx <- max(tt[, j]); rng <- mx - mn
    if (rng > 0) tt[, j] <- (tt[, j] - mn) / rng
  }
  tt
}

dist_rows <- function(t) {
  tt <- if (is.vector(t)) matrix(as.numeric(t), ncol = 1) else as.matrix(t)
  as.matrix(stats::dist(tt))
}

nloc <- function(t) if (is.matrix(t)) nrow(t) else length(t)

ensure_dyadic_J <- function(P, J) {
  if (is.null(J)) J <- log2(P)
  J_int <- as.integer(round(J))
  if (abs(P - 2^J_int) > .Machine$double.eps * max(1, P)) {
    stop(sprintf("P must be 2^J (dyadic). Got P=%s, J≈%.6f (rounded J=%d gives 2^J=%s).",
                 P, J, J_int, 2^J_int))
  }
  J_int
}

logit <- function(x) log(x/(1-x))
invlogit <- function(z) 1/(1+exp(-z))
`%||%` <- function(a, b) if (is.null(a)) b else a

pack_L <- function(L) L[lower.tri(L, diag = TRUE)]
unpack_L <- function(theta, m) {
  L <- matrix(0, m, m)
  L[lower.tri(L, diag = TRUE)] <- theta
  diag(L) <- abs(diag(L)) + 1e-8
  L
}

as_num_mat <- function(A) {
  if (!is.matrix(A)) A <- as.matrix(A)
  storage.mode(A) <- "double"
  A
}

normalize_t <- function(t, P) {
  if (is.null(t)) stop("t is NULL")
  if (is.vector(t)) {
    if (length(t) != P) stop(sprintf("t vector has length %d but P=%d", length(t), P))
    return(as.numeric(t))
  }
  tt <- as.matrix(t)
  if (nrow(tt) == P) return(tt)
  if (ncol(tt) == P) return(t(tt))
  stop(sprintf("t has incompatible shape: %d x %d; need P x d or length-P (P=%d).",
               nrow(tt), ncol(tt), P))
}

# Stick-breaking ----------------------------------------------------------
stick_to_pi <- function(v) {
  K <- length(v); pi <- numeric(K); cum <- 1
  for (k in 1:K) { pi[k] <- v[k] * cum; cum <- cum * (1 - v[k]) }
  pi
}
extend_sticks_until <- function(v, alpha, threshold) {
  tail <- prod(1 - v)
  while (tail > threshold) { v_new <- rbeta(1, 1, alpha); v <- c(v, v_new); tail <- tail * (1 - v_new) }
  v
}
update_v_given_z <- function(v, z, alpha) {
  K <- length(v); n_k <- tabulate(z, nbins = K); n_tail <- rev(cumsum(rev(n_k)))
  for (k in 1:K) { a <- 1 + n_k[k]; b <- alpha + if (k < K) n_tail[k+1] else 0; v[k] <- rbeta(1, a, b) }
  v
}

# Eigen helpers -----------------------------------------------------------
eig_Kx <- function(Kx) {
  ee <- eigen(Kx, symmetric = TRUE, only.values = FALSE)
  list(V = ee$vectors, s = pmax(ee$values, 0))
}
eig_Bshape <- function(L, M) {
  Bshape <- tcrossprod(L)
  trB <- sum(diag(Bshape))
  if (trB > 0) Bshape <- Bshape * (M / trB) # unit-trace
  ee <- eigen(Bshape, symmetric = TRUE)
  list(U = ee$vectors, lam = pmax(ee$values, 0))
}
project_curve <- function(Yi, mu, V) {
  Yi <- as_num_mat(Yi); mu <- as_num_mat(mu)
  if (is.null(V)) stop("Eigenvectors V are NULL.")
  crossprod(V, Yi - mu)
}
