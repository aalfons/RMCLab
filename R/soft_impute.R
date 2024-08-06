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

soft_impute <- function(X, lambda = fraction_grid(reverse = TRUE), 
                        relative = TRUE, rank.max = NULL, 
                        type = c("svd", "als"), thresh = 1e-05, 
                        maxit = 100L, trace.it = FALSE, 
                        final.svd = TRUE, 
                        # discretization of the imputed matrix is only done for 
                        # fitting the algorithm for the optimal lambda after 
                        # tuning, but it is currently ignored if a vector of 
                        # lambda is supplied
                        discretize = FALSE, values = NULL) {
  
  # initializations
  X <- as.matrix(X)
  
  # check arguments
  lambda <- sort(unique(lambda), decreasing = TRUE)
  nb_lambda <- length(lambda)
  relative <- isTRUE(relative)
  type <- match.arg(type)
  if (is.null(rank.max)) {
    # For the SVD algorithm, don't apply a rank restriction by default so that 
    # the solution solves the nuclear-norm convex matrix-completion problem. 
    # For the ALS algorithm, use the default of softImpute() (which is 2), as
    # this algorithm anyway only gives a local minimum of the nuclear-norm 
    # convex matrix-completion problem.
    rank.max <- if (type == "svd") min(dim(X)) - 1L else 2L
  }
  
  # center the data matrix
  X_centered <- softImpute::biScale(X, 
                                    row.center = FALSE, 
                                    row.scale = FALSE, 
                                    col.center = TRUE, 
                                    col.scale = FALSE)
  
  # if relative grid of tuning parameter values is requested, estimate the 
  # smallest value that sets all imputed values to 0 in the centered matrix
  if (relative) {
    d_max <- softImpute::lambda0(X_centered)  # largest singular value
  } else d_max <- 1
  
  # apply softImpute
  if (nb_lambda == 1L) {
    # apply softImpute with centered incomplete data matrix
    fit <- softImpute::softImpute(X_centered, rank.max = rank.max, 
                                  lambda = lambda * d_max, type = type, 
                                  thresh = thresh, maxit = maxit, 
                                  trace.it = trace.it, 
                                  final.svd = final.svd)
    # note that function complete() requires to use the data matrix on the 
    # original scale, only the predicted values for the missing elements are 
    # transformed back with argument 'unscale'
    X_imputed <- softImpute::complete(X, fit, unscale = TRUE)
  } else {
    # loop over values of the regularization parameter
    warm_start <- NULL
    X_imputed <- fit <- vector("list", length = nb_lambda)
    for (l in seq_along(lambda)) {
      # apply softImpute with centered incomplete data matrix
      fit[[l]] <- softImpute::softImpute(X_centered, rank.max = rank.max, 
                                         lambda = lambda[l] * d_max, 
                                         type = type, thresh = thresh, 
                                         maxit = maxit, trace.it = trace.it, 
                                         warm.start = warm_start, 
                                         final.svd = final.svd)
      # note that function complete() requires to use the data matrix on the 
      # original scale, only the predicted values for the missing elements are 
      # transformed back with argument 'unscale'
      X_imputed[[l]] <- softImpute::complete(X, fit[[l]], unscale = TRUE)
      # update starting values for next iteration
      warm_start <- fit[[l]]
    }
  }
  
  # construct object to be returned
  out <- list(lambda = lambda, lambda0 = d_max, svd = fit, X = X_imputed)
  
  # if requested, add discretized imputed matrix 
  # (only for fitting the algorithm for the optimal lambda after tuning)
  if (nb_lambda == 1L && isTRUE(discretize)) {
    # if not supplied, obtain values of rating scale (categories)
    if (is.null(values)) values <- unique(X[!is.na(X)])
    # obtain minimum and maximum value for discretization
    min_value <- min(values)
    max_value <- max(values)
    # discretize the imputed matrix
    # FIXME: this only works for integer sequence from min_value to max_value
    out$X_discretized <- pmin(pmax(round(X_imputed), min_value), max_value)
  }
  
  # add class and return object
  class(out) <- "soft_impute"
  out
}
