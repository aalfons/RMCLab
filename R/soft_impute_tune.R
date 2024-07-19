# *************************************
# Authors: Andreas Alfons
#          Erasmus University Rotterdam
# *************************************


## control object for selection of tuning parameter as a fraction of lambda0()
#' @export

fraction_control <- function(start = 0.01, end = 1, nb_lambda = 10L) {
  if (start <= 0) {
    stop("starting fraction of the tuning parameter must be greater than 0")
  }
  if (end > 1) {
    stop("ending fraction of the tuning parameter must be smaller than or equal to 1")
  }
  if (nb_lambda <= 0) {
    stop("number of tuning parameter values must be greater than 0")
  }
  out <- list(start = start, end = end, nb_lambda = nb_lambda)
  class(out) <- "fraction_control"
  out
}


## obtain grid of tuning parameter values

get_grid <- function(object, ...) UseMethod("get_grid")

get_grid.fraction_control <- function(object, ...) {
  seq_log <- seq(from = log(object$start), to = log(object$end), 
                 length.out = object$nb_lambda)
  exp(seq_log)
}


## function for tuning the penalty parameter via data splitting strategies
#' @export

soft_impute_tune <- function(X, lambda = fraction_control(), 
                             splits = holdout_control(), ..., 
                             discretize = TRUE, values = NULL) {
  
  # initializations
  X <- as.matrix(X)
  
  # check values of tuning parameter
  if (inherits(lambda, "fraction_control")) {
    # obtain grid of fractions for the tuning parameter
    pct_lambda <- get_grid(lambda)
    # center the data matrix
    X_centered <- softImpute::biScale(X, 
                                      row.center = FALSE, 
                                      row.scale = FALSE, 
                                      col.center = TRUE, 
                                      col.scale = FALSE)
    # estimate the smallest lambda that sets everything to zero and 
    # obtain the final grid of tuning parameter values
    lambda <- pct_lambda * softImpute::lambda0(X_centered)
  } else {
    # ensure values of tuning parameter are sorted
    lambda <- sort(unique(lambda))
  }
  if (length(lambda) == 1L) {
    stop("only one value of 'lambda'; use function soft_impute() instead")
  }
  
  # construct index vector of observed values
  observed <- which(!is.na(X))
  
  # create splits for tuning parameter validation
  if (inherits(splits, "split_control")) {
    splits <- create_splits(observed, control = splits)
  }
  
  # apply soft_impute() to the different training data sets
  fit_train <- lapply(splits, function(indices, ...) {
    # create training data where elements from test set are set to NA
    X_train <- X
    X_train[indices] <- NA_real_
    # apply soft_impute() to training data
    fit_train <- soft_impute(X_train, lambda = lambda, ...)
  }, ...)

  # extract predictions for the elements in the different test sets and compute 
  # prediction loss
  tuning_loss <- mapply(function(indices, fit) {
    # extract elements from test set as a vector
    X_test <- X[indices]
    # compute prediction loss (squared error)
    sapply(fit$X, function(X_hat) (X_test - X_hat[indices])^2)
  }, indices = splits, fit = fit_train, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  # combine results for tuning parameter validation into one matrix
  # (each column holds the prediction losses for the corresponding lambda)
  tuning_loss <- do.call(rbind, tuning_loss)
  # compute column means (mean squared error)
  tuning_loss <- colMeans(tuning_loss)
  
  # select the optimal lambda: revert the vectors so that in the unlikely 
  # case of ties, we select the lambda with stronger penalization
  which_opt <- rev(seq_along(lambda))[which.min(rev(tuning_loss))]
  lambda_opt <- lambda[which_opt]
  
  # apply soft_impute() with optimal tuning parameter
  # TODO: this could be made more efficient by adding an argument 'centered' 
  #       to soft_impute() in case the data matrix is already centered
  fit_opt <- soft_impute(X, lambda = lambda_opt, ..., discretize = discretize, 
                         values = values)

  ## construct list of relevant output
  out <- list(lambda = lambda, tuning_loss = tuning_loss, 
              lambda_opt = lambda_opt, fit = fit_opt)
  class(out) <- "soft_impute_tuned"
  out
  
}
