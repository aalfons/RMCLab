# **************************************
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# **************************************

## function for tuning the penalty parameter via data splitting strategies
#' @export

rdmc_tune <- function(X, values, lambda, 
                    loss = c("pseudo_huber", "absolute", "bounded"), 
                    loss_tuning = NULL, splits = holdout_control(), ...) {
  
  # initializations
  X <- as.matrix(X)
  
  # check arguments
  values <- sort(values)          # ensure values of rating scale are sorted
  lambda <- sort(unique(lambda))  # ensure values of tuning parameter are sorted
  loss <- match.arg(loss)
  
  # check bound in case bounded loss function
  if (is.null(loss_tuning)) {
    loss_tuning <- switch(loss, 
                          pseudo_huber = 1,
                          absolute = NA_real_,
                          bounded = (max(values) - min(values)) / 2)
  }
  
  if (length(lambda) == 1L) {
    stop("only one value of 'lambda'; use function rdmc() instead")
  }
  
  # create splits for tuning parameter validation
  observed <- which(!is.na(X))  # returns vector, which is what we need
  if (inherits(splits, "split_control")) {
    splits <- create_splits(observed, control = splits)
  }

  # # perform tuning parameter validation
  # tuning_loss <- lapply(splits, function(indices, ...) {
  #   # extract elements from test set as a vector
  #   X_test <- X[indices]
  #   # create training data where elements from test set are set to NA
  #   X_train <- X
  #   X_train[indices] <- NA_real_
  #   # apply robust discrete matrix completion to training data
  #   fit_train <- rdmc(X_train, values = values, lambda = lambda, loss = loss,
  #                     loss_tuning = loss_tuning, ...)
  #   # extract predictions for the elements in the test set and compute 
  #   # prediction loss
  #   sapply(fit_train$X, function(X_hat) {
  #     if (loss == "pseudo_huber") {
  #       pseudo_huber(X_test - X_hat[indices], delta = loss_tuning)
  #     } else if (loss == "absolute") {
  #       abs(X_test - X_hat[indices])
  #     } else {
  #       bounded(X_test - X_hat[indices], bound = loss_tuning)
  #     }
  #   })
  # })
  
  
  # fit robust discrete matrix completion to the different training data sets
  fit_train <- lapply(splits, function(indices, ...) {
    # create training data where elements from test set are set to NA
    X_train <- X
    X_train[indices] <- NA_real_
    # apply robust discrete matrix completion to training data
    fit_train <- rdmc(X_train, values = values, lambda = lambda, loss = loss,
                      loss_tuning = loss_tuning, ...)
  })
  
  # extract predictions for the elements in the different test sets and compute 
  # prediction loss
  tuning_loss <- mapply(function(indices, fit) {
    # extract elements from test set as a vector
    X_test <- X[indices]
    # compute prediction loss
    sapply(fit$X, function(X_hat) {
      if (loss == "pseudo_huber") {
        pseudo_huber(X_test - X_hat[indices], delta = loss_tuning)
      } else if (loss == "absolute") {
        abs(X_test - X_hat[indices])
      } else {
        bounded(X_test - X_hat[indices], bound = loss_tuning)
      }
    })
    
  }, indices = splits, fit = fit_train, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  # combine results for tuning parameter validation into one data frame
  # (each column holds the prediction losses for the corresponding lambda)
  tuning_loss <- do.call(rbind, tuning_loss)
  # compute column means
  tuning_loss <- colMeans(tuning_loss)
  
  # select the optimal lambda: revert the vectors so that in the unlikely 
  # case of ties, we select the lambda with stronger penalization
  which_opt <- rev(seq_along(lambda))[which.min(rev(tuning_loss))]
  lambda_opt <- lambda[which_opt]
  
  # compute starting values for L, Z, and Theta by aggregating solutions 
  # obtained from the different training data
  L <- sapply(fit_train, function(fit) fit$L[[which_opt]], 
              simplify = "array")
  L <- apply(L, 1:2, local_mode)
  Z <- sapply(fit_train, function(fit) fit$Z[[which_opt]], 
              simplify = "array")
  Z <- apply(Z, 1:2, mean)
  Theta <- sapply(fit_train, function(fit) fit$Theta[[which_opt]], 
                  simplify = "array")
  Theta <- apply(Theta, 1:2, mean)
  
  # apply robust discrete matrix completion with optimal tuning parameter
  fit_opt <- rdmc(X, values = values, lambda = lambda_opt, loss = loss,
                  loss_tuning = loss_tuning, ..., L = L, Z = Z, Theta = Theta)

  # construct list of relevant output
  out <- list(lambda = lambda, tuning_loss = tuning_loss, 
              lambda_opt = lambda_opt, final = fit_opt)
  class(out) <- "rdmc_tuned"
  out
  
}


# function to compute mode (for starting value of L in fit with optimal lambda)
local_mode <- function(x) {
  out <- DescTools::Mode(x)  # TODO: avoid this dependency
  if (length(out) > 1) sample(out, 1)
  else out
}
