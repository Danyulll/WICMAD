
# Utilities ---------------------------------------------------------------

#' Scale Input to [0,1] Range
#'
#' Scales input coordinates to the unit interval \[0,1\] for each dimension.
#'
#' @param t Numeric vector or matrix of input coordinates
#' @return Scaled coordinates in \[0,1\] range
#' @examples
#' scale_t01(c(1, 3, 5, 7))  # Returns c(0, 0.33, 0.67, 1)
#' scale_t01(matrix(c(1,2,3,4), ncol=2))  # Scales each column separately
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

#' Compute Distance Matrix Between Rows
#'
#' Computes pairwise distances between rows of input coordinates.
#'
#' @param t Numeric vector or matrix of input coordinates
#' @return Distance matrix between all pairs of rows
#' @examples
#' dist_rows(c(1, 2, 3))  # 3x3 distance matrix
dist_rows <- function(t) {
  tt <- if (is.vector(t)) matrix(as.numeric(t), ncol = 1) else as.matrix(t)
  as.matrix(stats::dist(tt))
}

#' Get Number of Locations
#'
#' Returns the number of input locations (rows for matrix, length for vector).
#'
#' @param t Numeric vector or matrix of input coordinates
#' @return Integer number of locations
#' @examples
#' nloc(c(1, 2, 3))  # Returns 3
#' nloc(matrix(1:6, ncol=2))  # Returns 3
nloc <- function(t) if (is.matrix(t)) nrow(t) else length(t)

#' Ensure Dyadic Wavelet Decomposition Level
#'
#' Ensures the number of wavelet decomposition levels J is appropriate for signal length P.
#'
#' @param P Integer signal length
#' @param J Integer or NULL number of decomposition levels
#' @return Integer number of decomposition levels
#' @examples
#' ensure_dyadic_J(64, NULL)  # Returns 6 (\eqn{\log_2(64)})
#' ensure_dyadic_J(64, 4)     # Returns 4
ensure_dyadic_J <- function(P, J) {
  if (is.null(J)) J <- log2(P)
  J_int <- as.integer(round(J))
  if (abs(P - 2^J_int) > .Machine$double.eps * max(1, P)) {
    stop(sprintf("P must be 2^J (dyadic). Got P=%s, J≈%.6f (rounded J=%d gives 2^J=%s).",
                 P, J, J_int, 2^J_int))
  }
  J_int
}

#' Logit Transformation
#'
#' Computes the logit transformation: \eqn{\text{logit}(x) = \log(x/(1-x))}
#'
#' @param x Numeric vector in (0,1)
#' @return Numeric vector of logit values
#' @examples
#' logit(0.5)  # Returns 0
#' logit(c(0.1, 0.9))  # Returns c(-2.2, 2.2)
logit <- function(x) log(x/(1-x))
#' Inverse Logit Transformation
#'
#' Computes the inverse logit transformation: \eqn{\text{invlogit}(z) = 1/(1+e^{-z})}
#'
#' @param z Numeric vector
#' @return Numeric vector in (0,1)
#' @examples
#' invlogit(0)  # Returns 0.5
#' invlogit(c(-2, 2))  # Returns c(0.12, 0.88)
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
  if (!is.matrix(A)) {
    A <- as.matrix(A)
  } else if (storage.mode(A) != "double") {
    storage.mode(A) <- "double"
  }
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
#' Stick-Breaking to Probability Vector
#'
#' Converts stick-breaking weights to probability vector using the stick-breaking process:
#' \eqn{\pi_k = v_k \prod_{j=1}^{k-1} (1-v_j)}
#'
#' @param v Numeric vector of stick-breaking weights in [0,1]
#' @return Numeric vector of probabilities that sum to 1
#' @examples
#' stick_to_pi(c(0.5, 0.3, 0.8))  # Returns c(0.5, 0.15, 0.12)
stick_to_pi <- function(v) {
  K <- length(v); pi <- numeric(K); cum <- 1
  for (k in 1:K) { pi[k] <- v[k] * cum; cum <- cum * (1 - v[k]) }
  pi
}
#' Extend Stick-Breaking Until Threshold
#'
#' Extends stick-breaking weights until the remaining mass falls below threshold.
#' Uses \eqn{v_k \sim \text{Beta}(1, \alpha)} for new sticks.
#'
#' @param v Numeric vector of existing stick weights
#' @param alpha Numeric concentration parameter for Beta distribution
#' @param threshold Numeric threshold for remaining mass
#' @return Extended vector of stick weights
extend_sticks_until <- function(v, alpha, threshold) {
  tail <- prod(1 - v)
  while (tail > threshold) { v_new <- rbeta(1, 1, alpha); v <- c(v, v_new); tail <- tail * (1 - v_new) }
  v
}
#' Update Stick Weights Given Cluster Assignments
#'
#' Updates stick-breaking weights given cluster assignments using Beta posteriors:
#' \eqn{v_k \sim \text{Beta}(1 + n_k, \alpha + \sum_{j=k+1}^K n_j)}
#'
#' @param v Numeric vector of stick weights
#' @param z Integer vector of cluster assignments
#' @param alpha Numeric concentration parameter
#' @return Updated vector of stick weights
update_v_given_z <- function(v, z, alpha) {
  K <- length(v); n_k <- tabulate(z, nbins = K); n_tail <- rev(cumsum(rev(n_k)))
  for (k in 1:K) { a <- 1 + n_k[k]; b <- alpha + if (k < K) n_tail[k+1] else 0; v[k] <- rbeta(1, a, b) }
  v
}

# Eigen helpers -----------------------------------------------------------
#' Eigendecomposition of Kernel Matrix
#'
#' Computes eigendecomposition of kernel matrix \eqn{K_x = U \Lambda U^T}.
#'
#' @param Kx Numeric kernel matrix
#' @return List with eigenvectors U and eigenvalues Lambda
eig_Kx <- function(Kx) {
  ee <- eigen(Kx, symmetric = TRUE, only.values = FALSE)
  list(V = ee$vectors, s = pmax(ee$values, 0))
}
#' Eigendecomposition of Coregionalization Matrix
#'
#' Computes eigendecomposition of coregionalization matrix \eqn{B = L L^T = U \Lambda U^T}.
#'
#' @param L Numeric lower triangular matrix
#' @param M Integer number of channels
#' @return List with eigenvectors U and eigenvalues Lambda
eig_Bshape <- function(L, M) {
  Bshape <- tcrossprod(L)
  trB <- sum(diag(Bshape))
  if (trB > 0) Bshape <- Bshape * (M / trB) # unit-trace
  ee <- eigen(Bshape, symmetric = TRUE)
  list(U = ee$vectors, lam = pmax(ee$values, 0))
}
#' Project Curve onto Subspace
#'
#' Projects curve onto subspace defined by eigenvectors: \eqn{Y_{\text{proj}} = V V^T (Y - \mu) + \mu}
#'
#' @param Yi Numeric matrix of curve data
#' @param mu Numeric vector of mean
#' @param V Numeric matrix of eigenvectors
#' @return Numeric matrix of projected curve
project_curve <- function(Yi, mu, V) {
  Yi <- as_num_mat(Yi); mu <- as_num_mat(mu)
  if (is.null(V)) stop("Eigenvectors V are NULL.")
  crossprod(V, Yi - mu)
}
