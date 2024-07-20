# ************************************
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ************************************


## function for tuning the penalty parameter via data splitting strategies
#' @export

rdmc_tune <- function(X, values = NULL, lambda = rel_lambda_control(), 
                      splits = holdout_control(), 
                      loss = c("bounded", "absolute", "pseudo_huber"),
                      loss_const = NULL, mu = 0.1, ...) {
  
  # initializations
  X <- as.matrix(X)
  
  # construct index vector of observed values
  observed <- which(!is.na(X))
  
  # check values of rating scale
  if (is.null(values)) values <- unique(X[observed])
  values <- sort(values)  # ensure values of rating scale are sorted
  # check values of tuning parameter
  if (inherits(lambda, "rel_lambda_control")) {
    # obtain grid of relative values for the tuning parameter
    rel_lambda <- get_grid(lambda)
    # center data matrix with columnwise median of observed data
    medians <- apply(X, 2L, median, na.rm = TRUE)
    X_centered <- sweep(X, 2, medians, FUN = "-")
    # obtain the final grid of tuning parameter values relative to the largest 
    # singular value of the centered data matrix with missing values replaced 
    # by 0 (note that in this case, this is not an estimate of the smallest 
    # lambda that sets everything to 0, but we use it here simply to relate 
    # the grid of lambda values to the specific data at hand)
    lambda <- rel_lambda * mu * softImpute::lambda0(X_centered)
  } else {
    # ensure values of tuning parameter are sorted
    lambda <- sort(unique(lambda), decreasing = FALSE)
  }
  if (length(lambda) == 1L) {
    stop("only one value of 'lambda'; use function rdmc() instead")
  }
  # have_autotune <- inherits(lambda, "autotune_control")
  # if (!have_autotune) {
  #   # ensure values of tuning parameter are sorted
  #   lambda <- sort(unique(lambda), decreasing = TRUE)
  #   if (length(lambda) == 1L) {
  #     stop("only one value of 'lambda'; use function rdmc() instead")
  #   }
  # }
  # check loss function
  loss <- match.arg(loss)
  if (is.null(loss_const)) {
    # set default constant for loss function (if applicable)
    loss_const <- switch(loss, bounded = (max(values) - min(values)) / 2, 
                         absolute = NA_real_, pseudo_huber = 1)
  }
  
  # create splits for tuning parameter validation
  if (inherits(splits, "split_control")) {
    splits <- create_splits(observed, control = splits)
  }
  
  # # if requested, automatically determine values of tuning parameter
  # if (have_autotune) {
  #   # iteratively fit robust discrete matrix completion with increasing values
  #   # until tuning parameter is large enough so that first soft-thresholded SVD 
  #   # step results in all singular values being zero
  #   fit <- rdmc(X, values = values, lambda = lambda, loss = loss,
  #               loss_const = loss_const, ...)
  #   # extract values of tuning parameter
  #   lambda <- fit$lambda
  # }

  # fit robust discrete matrix completion to the different training data sets
  fit_train <- lapply(splits, function(indices, ...) {
    # create training data where elements from test set are set to NA
    X_train <- X
    X_train[indices] <- NA_real_
    # apply robust discrete matrix completion to training data
    fit_train <- rdmc(X_train, values = values, lambda = lambda, loss = loss,
                      loss_const = loss_const, mu = mu, ...)
  }, ...)
  
  # extract predictions for the elements in the different test sets and compute 
  # prediction loss
  tuning_loss <- mapply(function(indices, fit) {
    # extract elements from test set as a vector
    X_test <- X[indices]
    # compute prediction loss
    sapply(fit$X, function(X_hat) {
      if (loss == "pseudo_huber") {
        pseudo_huber(X_test - X_hat[indices], delta = loss_const)
      } else if (loss == "absolute") {
        abs(X_test - X_hat[indices])
      } else {
        bounded(X_test - X_hat[indices], bound = loss_const)
      }
    })
    
  }, indices = splits, fit = fit_train, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  # combine results for tuning parameter validation into one matrix
  # (each column holds the prediction losses for the corresponding lambda)
  tuning_loss <- do.call(rbind, tuning_loss)
  # compute column means
  tuning_loss <- colMeans(tuning_loss)
  
  # select the optimal lambda: reverse vectors so that in the unlikely case of 
  # ties, we select the lambda with stronger penalization
  which_opt <- rev(seq_along(lambda))[which.min(rev(tuning_loss))]
  lambda_opt <- lambda[which_opt]
  
  # # prepare output
  # if (have_autotune) {
  #   
  #   # construct list of relevant output
  #   out <- list(lambda = lambda, tuning_loss = tuning_loss, 
  #               which_opt = which_opt, lambda_opt = lambda_opt, 
  #               fit = fit)
  #   class(out) <- c("rdmc_autotuned", "rdmc_tuned")
  #   
  # } else {
    
    # Note: It's possible that on different training sets, we get different 
    # medians in some variables so that the discrete constraint is different 
    # for those variables. But for the starting values, it shouldn't matter 
    # (L doesn't even have to satisfy any discrete constraint), as we simply 
    # use L and Theta for the first soft-thresholding SVD step to update Z, 
    # and subsequently L is updated to satisfy the discrete constrains using 
    # the medians on the full data set. Hence we could just take the mean or 
    # median over the training sets for L rather than the mode (which currently 
    # uses random sampling in case of multiple modes). Specifically, we apply 
    # the soft-thresholding SVD step to L + Theta/mu, so it may be smoother to 
    # use the linearity of the mean to compute both starting values for L and 
    # Theta, rather than destroying the linear relationship by using the mode 
    # for L and the mean for Theta.
    
    # compute starting values for L and Theta by aggregating solutions obtained 
    # from the different training data
    L <- sapply(fit_train, function(fit) fit$L[[which_opt]], 
                simplify = "array")
    L <- apply(L, 1:2, median)
    Theta <- sapply(fit_train, function(fit) fit$Theta[[which_opt]], 
                    simplify = "array")
    Theta <- apply(Theta, 1:2, mean)
    
    # apply robust discrete matrix completion with optimal tuning parameter
    # TODO: this could be made more efficient by adding arguments to rdmc() 
    #       in case the data matrix is already centered
    fit_opt <- rdmc(X, values = values, lambda = lambda_opt, loss = loss,
                    loss_const = loss_const, mu = mu, ..., L = L, 
                    Theta = Theta)
    
    # construct list of relevant output
    out <- list(lambda = lambda, tuning_loss = tuning_loss, 
                lambda_opt = lambda_opt, fit = fit_opt)
    class(out) <- "rdmc_tuned"
  
  # }
  
  # return output
  out
  
}


# ## function to compute mode for starting value of L in fit with optimal lambda
# local_mode <- function(x) {
#   # compute contingency table
#   tab <- table(x)
#   # determine the maximum frequency
#   max <- max(tab)
#   # keep all values that occur with the maximum frequency
#   out <- as.numeric(names(tab)[tab == max])
#   # if there are multiple modes, randomly select one to be returned
#   if (length(out) > 1) sample(out, 1)
#   else out
# }
