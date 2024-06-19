# *************************************
# Authors: Andreas Alfons
#          Erasmus University Rotterdam
#
#          Aurore Archimbaud
#          Toulouse Business School
# *************************************


## wrapper function for softImpute() that returns continuous and discretized
## predictions for a vector of values for the regularization parameter
#' @importFrom softImpute biScale complete softImpute
#' @export
soft_impute <- function(X, lambda, rank.max = 2, type = "svd", thresh = 1e-05, 
                        maxit = 100, trace.it = FALSE, final.svd = TRUE, 
                        # discretization of the imputed matrix is only done for 
                        # fitting the algorithm for the optimal lambda after 
                        # tuning, but it is currently ignored if a vector of 
                        # lambda is supplied
                        discretize = FALSE, values = NULL) {
  
  # initializations
  X <- as.matrix(X)
  
  # check arguments
  lambda <- sort(unique(lambda)) # ensure values of tuning parameter are sorted
  n_lambda <- length(lambda)
  
  # center the data matrix
  X_centered <- softImpute::biScale(X, 
                                    row.center = FALSE, 
                                    row.scale = FALSE, 
                                    col.center = TRUE, 
                                    col.scale = FALSE)
  
  # apply softImpute
  if (n_lambda == 1L) {
    # apply softImpute with centered incomplete data matrix
    fit <- softImpute::softImpute(X_centered, rank.max = rank.max, 
                                  lambda = lambda, type = type, 
                                  thresh = thresh, maxit = maxit, 
                                  trace.it = trace.it, 
                                  final.svd = final.svd)
    # note that function complete() requires to use the data matrix on the 
    # original scale, only the predicted values for the missing cells are 
    # transformed back with argument 'unscale'
    X_imputed <- softImpute::complete(X, fit, unscale = TRUE)
  } else {
    # loop over values of the regularization parameter
    warm_start <- NULL
    X_imputed <- fit <- vector("list", length = length(lambda))
    for (l in seq_along(lambda)) {
      # apply softImpute with centered incomplete data matrix
      fit[[l]] <- softImpute::softImpute(X_centered, rank.max = rank.max, 
                                         lambda = lambda[l], type = type, 
                                         thresh = thresh, maxit = maxit,
                                         trace.it = trace.it, 
                                         warm.start = warm_start, 
                                         final.svd = final.svd)
      # note that function complete() requires to use the data matrix on the 
      # original scale, only the predicted values for the missing cells are 
      # transformed back with argument 'unscale'
      X_imputed[[l]] <- softImpute::complete(X, fit[[l]], unscale = TRUE)
      # update starting values for next iteration
      warm_start <- fit[[l]]
    }
  }
  
  # construct object to be returned
  out <- list(lambda = lambda, svd = fit, X = X_imputed)
  
  # if requested, add discretized imputed matrix 
  # (only for fitting the algorithm for the optimal lambda after tuning)
  if (n_lambda == 1L && isTRUE(discretize)) {
    # if not supplied, obtain values of rating scale (categories)
    if (is.null(values)) values <- unique(X[!is.na(X)])
    # obtain minimum and maximum value for discretization
    min_value <- min(values)
    max_value <- max(values)
    # discretize the imputed matrix
    out$X_discretized <- pmin(pmax(round(X_imputed), min_value), max_value)
  }
  
  # add class and return object
  class(out) <- "soft_impute"
  out
}
