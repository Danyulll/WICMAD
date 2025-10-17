
# Wavelet forward/inverse & precomputation -------------------------------

#' Forward 1D Wavelet Transform
#'
#' Performs forward discrete wavelet transform on a 1D signal using the waveslim package.
#' Returns wavelet coefficients in a flattened vector format with indexing information.
#'
#' @param y Numeric vector of length P representing the input signal
#' @param wf Character string specifying the wavelet filter (default: "la8")
#' @param J Integer specifying the number of decomposition levels (default: NULL, auto-determined)
#' @param boundary Character string specifying boundary condition (default: "periodic")
#'
#' @return List containing:
#'   \item{coeff}{Numeric vector of flattened wavelet coefficients}
#'   \item{map}{List containing transform metadata (J, wf, boundary, P, idx)}
#'
#' @details The function uses the waveslim package to perform the wavelet transform
#' and reorganizes the output into a flattened coefficient vector with proper indexing.
#' The map object contains all necessary information for reconstruction.
#'
#' @examples
#' y <- rnorm(64)
#' result <- wt_forward_1d(y, wf = "la8", J = 4)
#' reconstructed <- wt_inverse_1d(result$coeff, result$map)
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
#' Inverse 1D Wavelet Transform
#'
#' Reconstructs original signal from wavelet coefficients using the waveslim package.
#' This is the inverse operation of wt_forward_1d().
#'
#' @param coeff_vec Numeric vector of flattened wavelet coefficients
#' @param map List containing transform metadata from wt_forward_1d()
#'
#' @return Numeric vector of length P representing the reconstructed signal
#'
#' @details The function reconstructs the original signal by:
#' 1. Reorganizing coefficients back into the waveslim dwt format
#' 2. Applying the inverse discrete wavelet transform
#' 3. Returning the reconstructed signal
#'
#' @examples
#' y <- rnorm(64)
#' result <- wt_forward_1d(y, wf = "la8", J = 4)
#' reconstructed <- wt_inverse_1d(result$coeff, result$map)
#' all.equal(y, reconstructed)  # Should be TRUE
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
#' Forward Wavelet Transform for Matrix
#'
#' Performs forward wavelet transform on each column of a matrix independently.
#' Useful for multivariate functional data where each column represents a different channel.
#'
#' @param y_mat Matrix of size P x M where P is the number of time points and M is the number of channels
#' @param wf Character string specifying the wavelet filter (default: "la8")
#' @param J Integer specifying the number of decomposition levels (default: NULL, auto-determined)
#' @param boundary Character string specifying boundary condition (default: "periodic")
#'
#' @return List of length M, where each element is the result of wt_forward_1d() for the corresponding column
#'
#' @details This function applies wt_forward_1d() to each column of the input matrix,
#' making it suitable for multivariate functional data analysis.
#'
#' @examples
#' y_mat <- matrix(rnorm(128), nrow = 64, ncol = 2)
#' result <- wt_forward_mat(y_mat, wf = "la8", J = 4)
wt_forward_mat <- function(y_mat, wf="la8", J=NULL, boundary="periodic") {
  M <- ncol(y_mat); out <- vector("list", M)
  for (m in 1:M) out[[m]] <- wt_forward_1d(y_mat[,m], wf, J, boundary)
  out
}
#' Precompute Wavelet Transforms
#'
#' Precomputes wavelet transforms for all curves in a list to avoid repeated computation during MCMC.
#'
#' @param Y_list List of matrices, where each matrix represents a functional curve
#' @param wf Character string specifying the wavelet filter
#' @param J Integer specifying the number of decomposition levels
#' @param boundary Character string specifying boundary condition
#'
#' @return List of the same length as Y_list, where each element contains the wavelet transform results
#'
#' @details This function is called once at the beginning of the MCMC algorithm to precompute
#' all wavelet transforms.
#'
#' @examples
#' Y_list <- list(matrix(rnorm(128), nrow = 64, ncol = 2), 
#'                matrix(rnorm(128), nrow = 64, ncol = 2))
#' precomp <- precompute_wavelets(Y_list, wf = "la8", J = 4, boundary = "periodic")
precompute_wavelets <- function(Y_list, wf, J, boundary) {
  lapply(Y_list, function(mat) wt_forward_mat(mat, wf, J, boundary))
}
#' Stack Wavelet Coefficients from Precomputed Transforms
#'
#' Extracts and stacks wavelet coefficients from precomputed transforms for a specific set of curves.
#' This function is used during MCMC to efficiently access wavelet coefficients for cluster updates.
#'
#' @param precomp List of precomputed wavelet transforms from precompute_wavelets()
#' @param idx Integer vector specifying which curves to extract (indices into precomp)
#' @param M Integer specifying the number of channels
#' @param bias_coeff Optional list of bias coefficients to subtract from each channel
#'
#' @return List containing:
#'   \item{D_arr}{3D array of size (ncoeff, N, M) containing stacked coefficients}
#'   \item{maps}{List of wavelet transform metadata for the first curve}
#'
#' @details This function efficiently extracts wavelet coefficients for a subset of curves
#' and organizes them into a 3D array structure suitable for MCMC updates. If bias coefficients
#' are provided, they are subtracted from the corresponding channels.
#'
#' @examples
#' # Assuming precomp is from precompute_wavelets()
#' result <- stack_D_from_precomp(precomp, idx = c(1, 3, 5), M = 2)
#' dim(result$D_arr)  # Should be (ncoeff, 3, 2)
stack_D_from_precomp <- function(precomp, idx, M, bias_coeff = NULL) {
  ncoeff <- length(precomp[[ idx[1] ]][[1]]$coeff)
  N <- length(idx)
  D_arr <- array(NA_real_, dim=c(ncoeff, N, M))
  for (jj in seq_along(idx)) {
    i <- idx[jj]
    for (m in 1:M) {
      D_arr[, jj, m] <- precomp[[i]][[m]]$coeff
      if (!is.null(bias_coeff)) {
        bc <- bias_coeff[[m]]
        if (!is.null(bc) && length(bc) == ncoeff) {
          D_arr[, jj, m] <- D_arr[, jj, m] - bc
        }
      }
    }
  }
  maps <- precomp[[ idx[1] ]]
  list(D_arr = D_arr, maps = maps)
}
#' Compute Mean Functions from Wavelet Coefficients
#'
#' Reconstructs mean functions from wavelet coefficients using inverse transform.
#' This function converts wavelet-domain parameters back to the time domain for visualization
#' and interpretation.
#'
#' @param beta_ch List of length M containing wavelet coefficients for each channel
#' @param wf Character string specifying the wavelet filter
#' @param J Integer specifying the number of decomposition levels
#' @param boundary Character string specifying boundary condition
#' @param P Integer specifying the number of time points in the original signal
#'
#' @return Matrix of size P x M containing the reconstructed mean functions for each channel
#'
#' @details This function is used to reconstruct mean functions from wavelet coefficients
#' stored during MCMC sampling. It creates a template using a zero matrix to ensure
#' proper wavelet structure, then applies the inverse transform to each channel's coefficients.
#'
#' @examples
#' # Assuming beta_ch contains wavelet coefficients for 2 channels
#' mu <- compute_mu_from_beta(beta_ch, wf = "la8", J = 4, boundary = "periodic", P = 64)
#' dim(mu)  # Should be (64, 2)
compute_mu_from_beta <- function(beta_ch, wf, J, boundary, P) {
  M <- length(beta_ch); zeros <- matrix(0, nrow=P, ncol=M)
  tmpl <- wt_forward_mat(zeros, wf, J, boundary); mu <- matrix(0, P, M)
  for (m in 1:M) mu[,m] <- wt_inverse_1d(beta_ch[[m]], tmpl[[m]]$map)
  mu
}
