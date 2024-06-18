# **************************************
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# **************************************

## control object for holdout validation
#' @export
holdout_control <- function(pct = 0.1, R = 10L) {
  if (pct <= 0 || pct >= 1) {
    stop("the percentage of holdout elements must be in the interval (0, 1)")
  }
  if (R < 1L) stop("holdout validation requires at least 1 replication")
  control <- list(pct = pct, R = R)
  class(control) <- "holdout_control"
  control
}

## control object for cross-validation
#' @export
fold_control <- function(K = 5L) {
  if (K < 2L) stop("cross-validation requires at least 2 folds")
  control <- list(K = K)
  class(control) <- "fold_control"
  control
}


## generic function for validation splits
#' @export
splits <- function(indices, control) UseMethod("splits", control)

## method to draw indices for holdout validation
#' @export
splits.holdout_control <- function(indices, control) {
  # initializations
  n <- length(indices)
  pct <- control$pct
  R <- control$R
  # check if enough indices have been supplied
  min <- ceiling(0.5000000001/pct) # avoid issues with floating point arithmetic
  if (n < min) {
    stop(sprintf("'indices' should contain at least %d elements", min))
  }
  # draw indices of holdout set
  m <- round(pct * n)
  replicate(R, sample(indices, size = m), simplify = FALSE)
}

## method to draw folds of indices for K-fold cross-validation
#' @export
splits.fold_control <- function(indices, control) {
  # initialization
  n <- length(indices)
  K <- control$K
  # check if enough indices have been supplied
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
  

## wrapper function to draw indices for holdout validation
#' @export
holdout <- function(indices, pct = 0.1, R = 10) {
  control <- holdout_control(pct = pct, R = R)
  splits(indices, control)
}


## wrapper function to draw folds of indices for K-fold cross-validation
#' @export
cv_folds <- function(indices, K = 5L) {
  control <- fold_control(K = K)
  splits(indices, control)
}


# ## function to draw indices for holdout validation
# holdout <- function(indices, pct = 0.1, R = 10) {
#   n <- length(indices)
#   min <- ceiling(0.5000000001/pct) # avoid issues with floating point arithmetic
#   if (n < min) {
#     stop(sprintf("'indices' should contain at least %d elements", min))
#   }
#   m <- round(pct * n)
#   replicate(R, sample(indices, size = m), simplify = FALSE)
# }
# 
# ## function to draw folds of indices for K-fold cross-validation
# cv_folds <- function(indices, K = 5L) {
#   n <- length(indices)
#   if (K < 2L) stop("cross-validation requires at least 2 folds")
#   if (n < K) stop(sprintf("'indices' should contain at least %d elements", K))
#   # permute observations
#   indices <- sample(indices)
#   # assign a block to each observation
#   blocks <- rep(seq_len(K), length.out = n)
#   # split the permuted observations according to the block they belong to
#   folds <- split(indices, blocks)
#   names(folds) <- NULL
#   folds
# }
