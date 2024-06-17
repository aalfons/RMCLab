# **************************************
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# **************************************

## function select the penalty parameter via K-fold cross-validation
#' @export

rdmc_cv <- function(X, values, lambda, 
                    loss = c("pseudo_huber", "absolute", "bounded"), 
                    loss_tuning = NULL, K = 5, folds = NULL, ...) {
  
  # check arguments
  values <- sort(values)       # make sure values of rating scale are sorted
  loss <- match.arg(loss)
  
  # check bound in case bounded loss function
  if (is.null(loss_tuning)) {
    loss_tuning <- switch(loss, 
                          pseudo_huber = 1,
                          absolute = NA_real_,
                          bounded = (max(values) - min(values)) / 2)
  }
  
  # create folds for K-fold cross-validation
  if (is.null(folds)) {
    n <- nrow(X)
    folds <- cv_folds(n, K = K)
  } else K <- length(folds)
  
  # perform cross-validation
  cv_results <- lapply(folds, function(indices, ...) {
    
    # split data matrix into training and test set
    X_train <- X[-indices, , drop = FALSE]
    X_test <- X[indices, , drop = FALSE]
    
    # apply robust discrete matrix completion to training data
    fit_train <- rdmc(X_train, values = values, lambda = lambda, loss = loss, 
                      loss_tuning = loss_tuning, ...)
    # FIXME: actually, we cannot obtain predictions for the removed rows
    # One solution would be not to use cross-validation, but repeatedly set 
    # some elements to be missing (perhaps with the restriction that at least 
    # one rating per user remains), and fit the algorithm to the nxp matrix 
    # that is a bit more incomplete.
  })
  
  
}
