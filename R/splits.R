# **************************************
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# **************************************

## function to draw indices for holdout validation
#' @export
holdout <- function(indices, pct = 0.1, R = 10) {
  n <- length(indices)
  min <- ceiling(0.5000000001/pct) # avoid issues with floating point arithmetic
  if (n < min) {
    stop(sprintf("'indices' should contain at least %d elements", min))
  }
  m <- round(pct * n)
  replicate(R, sample(indices, size = m), simplify = FALSE)
}

## function to draw folds of indices for K-fold cross-validation
#' @export
cv_folds <- function(indices, K = 5L) {
  n <- length(indices)
  if (K < 2L) stop("cross-validation requires at least 2 folds")
  if (n < K) stop(sprintf("'indices' should contain at least %d elements", K))
  # permute observations
  indices <- sample(indices)
  # assign a block to each observation
  blocks <- rep(seq_len(K), length.out = n)
  # split the permuted observations according to the block they belong to
  folds <- split(indices, blocks)
  names(folds) <- NULL
  folds
}
