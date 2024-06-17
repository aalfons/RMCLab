# **************************************
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# **************************************

#' @export

# Function to set up folds for K-fold cross-validation
cv_folds <- function(n, K = 5L) {
  # permute observations
  indices <- sample.int(n)
  # assign a block to each observation
  blocks <- rep(seq_len(K), length.out = n)
  # split the permuted observations according to the block they belong to
  folds <- split(indices, blocks)
  names(folds) <- NULL
  folds
}
