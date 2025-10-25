
# ICM likelihood via block cache -----------------------------------

#' Build ICM Cache for Efficient Likelihood Computation
#'
#' Builds or updates a cache for efficient ICM likelihood computation using eigendecomposition
#' and block diagonal structure. The cache stores precomputed eigendecompositions and
#' Cholesky factors to avoid repeated computation during MCMC.
#'
#' @param t Matrix of input locations
#' @param kern_cfg Kernel configuration list
#' @param kp Kernel parameters
#' @param L Lower triangular matrix for coregionalization
#' @param eta Vector of channel-specific noise parameters
#' @param tau_B Global scaling parameter
#' @param cache Existing cache object (optional, for updates)
#'
#' @return Updated cache object containing:
#'   \item{Ux_x}{Eigenvectors from kernel eigendecomposition}
#'   \item{lam_x}{Eigenvalues from kernel eigendecomposition}
#'   \item{Bshape}{Coregionalization matrix (unit trace)}
#'   \item{chol_list}{List of Cholesky factors for each eigenvalue block}
#'   \item{logdet_sum}{Sum of log determinants for likelihood computation}
#'   \item{key_kx}{Cache key for kernel eigendecomposition}
#'   \item{key_B}{Cache key for coregionalization matrix}
#'
#' @details The function implements: \eqn{\Sigma_j = \tau_B \lambda_j B + D_\eta},
#' where \eqn{\lambda_j} are eigenvalues from kernel eigendecomposition, \eqn{B} is the coregionalization
#' matrix, and \eqn{D_\eta} contains channel-specific noise. The cache is updated only when
#' necessary based on parameter changes.
#'
#' @examples
#' # Build initial cache
#' cache <- .build_icm_cache(t, kern_cfg, kp, L, eta, tau_B)
#' 
#' # Update cache with new parameters
#' cache <- .build_icm_cache(t, kern_cfg, new_kp, L, eta, tau_B, cache)
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

#' ICM Likelihood Computation for Single Curve
#'
#' Computes the log-likelihood for a single curve using precomputed cache from
#' .build_icm_cache(). Uses the block diagonal structure for efficient computation.
#'
#' @param y_resid Matrix of residuals (P x M) where P is number of time points and M is number of channels
#' @param cache Cache object from .build_icm_cache() containing precomputed eigendecompositions and Cholesky factors
#'
#' @return Numeric scalar representing the log-likelihood value
#'
#' @details The function computes the log-likelihood using the formula:
#' -0.5 * (P * M * log(2π) + logdet_sum + quadratic_form)
#' where the quadratic form is computed efficiently using precomputed Cholesky factors
#' for each eigenvalue block in the block diagonal structure.
#'
#' @examples
#' # Assuming cache is built and y_resid is available
#' loglik <- fast_icm_loglik_curve(y_resid, cache)
fast_icm_loglik_curve <- function(y_resid, cache) {
  if (is.null(cache$Ux_x) || is.null(cache$chol_list))
    stop("Cache not built: missing Ux_x or chol_list.")
  P_y  <- nrow(y_resid)
  P_ch <- length(cache$chol_list)
  if (P_y != P_ch)
    stop(sprintf("Cache dimension mismatch: y rows=%d, chol blocks=%d.", P_y, P_ch))

  # Use C++ implementation for maximum performance
  if (requireNamespace("RcppEigen", quietly = TRUE)) {
    tryCatch({
      return(fast_icm_loglik_curve_eigen(y_resid, cache$Ux_x, cache$chol_list, cache$logdet_sum))
    }, error = function(e) {
      warning("C++ implementation failed, falling back to R: ", e$message)
    })
  }
  
  # Fallback to optimized R implementation
  Ytil <- t(cache$Ux_x) %*% y_resid
  quad <- 0
  # Pre-allocate vectors to avoid repeated memory allocation
  M <- ncol(y_resid)
  v <- numeric(M)
  for (j in 1:P_ch) {
    Lj <- cache$chol_list[[j]]
    # Reuse pre-allocated vector instead of creating new one each iteration
    v[] <- Ytil[j, ]  # Fill existing vector instead of as.numeric(Ytil[j, ])
    w  <- backsolve(Lj, forwardsolve(t(Lj), v, upper.tri=TRUE, transpose=TRUE))
    quad <- quad + sum(w*w)
  }
  -0.5 * (P_ch * ncol(y_resid) * log(2*base::pi) + cache$logdet_sum + quad)
}
