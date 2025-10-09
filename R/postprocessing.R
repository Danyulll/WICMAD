
# Dahl / MAP helpers -----------------------------------------------------

dahl_partition <- function(Z) {
  if (is.null(dim(Z)) || length(dim(Z)) != 2) stop("Z must be an S x N matrix of labels.")
  S <- nrow(Z); N <- ncol(Z); PSM <- matrix(0, N, N)
  for (s in 1:S) { zs <- Z[s, ]; A <- outer(zs, zs, FUN = "==") * 1L; PSM <- PSM + A }
  PSM <- PSM / S; score <- numeric(S)
  for (s in 1:S) { zs <- Z[s, ]; A <- outer(zs, zs, FUN = "==") * 1L; score[s] <- sum((A - PSM)^2) }
  s_hat <- which.min(score); z_hat <- as.integer(factor(Z[s_hat, ])); K_hat <- length(unique(z_hat))
  list(z_hat = z_hat, K_hat = K_hat, s_hat = s_hat, PSM = PSM, score = score)
}
dahl_from_res <- function(res) { if (is.null(res$Z)) stop("res$Z not found"); dahl_partition(res$Z) }
.canon_labels <- function(z) { u <- unique(z); match(z, u) }
map_partition <- function(Z) {
  if (is.null(dim(Z)) || length(dim(Z)) != 2) stop("Z must be an S x N matrix of labels.")
  S <- nrow(Z); N <- ncol(Z); keys <- character(S)
  for (s in 1:S) keys[s] <- paste(.canon_labels(Z[s, ]), collapse = "-")
  tab <- table(keys); key_hat <- names(which.max(tab)); s_hat <- which(keys == key_hat)[1]
  z_hat <- .canon_labels(Z[s_hat, ]); K_hat <- length(unique(z_hat))
  list(z_hat = z_hat, K_hat = K_hat, s_hat = s_hat, key_hat = key_hat, freq = max(tab))
}
map_from_res <- function(res) { if (is.null(res$Z)) stop("res$Z not found"); map_partition(res$Z) }
