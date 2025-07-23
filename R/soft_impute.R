# *************************************
# Authors: Andreas Alfons
#          Erasmus University Rotterdam
#
#          Aurore Archimbaud
#          Toulouse Business School
# *************************************


## wrapper function for softImpute() that returns continuous and discretized
## predictions for a vector of values for the regularization parameter

#' Matrix completion via nuclear-norm regularization
#' 
#' Convenience wrapper for \code{\link[softImpute]{softImpute}()} that allows 
#' to supply a grid of values for the regularization parameter.  Other 
#' noteworthy differences with the original function are that the columns of 
#' the data matrix are centered internally, that some of the default values are 
#' different, and that the output is structured differently.  Moreover, in case 
#' of discrete rating-scale data, the wrapper function allows to include a 
#' discretization step after fitting the algorithm to map the imputed values to 
#' the rating scale of the observed values.
#' 
#' @param X  a matrix or data frame with missing values.
#' @param lambda  a numeric vector giving values of the regularization 
#' parameter.  See \code{\link{fraction_grid}()} for the default values.
#' @param relative  a logical indicating whether the values of the 
#' regularization parameter should be considered relative to a certain 
#' reference value computed from the data at hand.  If \code{TRUE} (the 
#' default), the values of \code{lambda} are multiplied with the value 
#' returned by \code{\link[softImpute]{lambda0}()}.
#' @param type  a character string specifying the type of algorithm. Possible 
#' values are \code{"svd"} and \code{"als"}. See 
#' \code{\link[softImpute]{softImpute}()} for details on the algorithms, but 
#' note that the default value here is \code{"svd"}.
#' @param rank.max  a positive integer giving a rank constraint.  See 
#' \code{\link[softImpute]{softImpute}()} for more details, but note that the 
#' default here is to use the minimum of the number of rows and columns minus 1 
#' if \code{type} is \code{"svd"}, and to use 2 if \code{type} is \code{"als"}.
#' @param thresh,maxit,trace.it,final.svd  see 
#' \code{\link[softImpute]{softImpute}()}.
#' @param discretize  a logical indicating whether to include a discretization 
#' step after fitting the algorithm (defaults to \code{FALSE}).  In case of 
#' discrete rating-scale data, this can be used to map the imputed values to 
#' the discrete rating scale of the observed values.
#' @param values  an optional numeric vector giving the possible values of 
#' discrete ratings.  This is ignored if \code{discretize} is \code{FALSE}.  
#' Currently, the possible values are assumed to be the same for all columns.  
#' If \code{NULL}, the unique values of the observed parts of \code{X} are 
#' used.
#' 
#' @return 
#' An object of class \code{"soft_impute"}.  The class structure is still 
#' experimental and may change.  
#' 
#' The following accessor functions are available:
#' \itemize{
#'   \item \code{\link{get_completed}()} to extract the completed (i.e., 
#'   imputed) data matrix with a specified value of the regularization 
#'   parameter,
#'   \item \code{\link{get_lambda}()} to extract the values of the 
#'   regularization parameter.
#' }
#' 
#' @author Andreas Alfons and Aurore Archimbaud
#' 
#' @references 
#' Hastie, T., Mazumder, R., Lee, J. D. and Zadeh, R. (2015) Matrix Completion 
#' and Low-Rank SVD via Fast Alternating Least Squares. \emph{Journal of 
#' Machine Learning Research}, \bold{16}(104), 3367--3402.
#' 
#' Mazumder, R., Hastie, T. and Tibshirani, R. (2010) Spectral Regularization 
#' Algorithms for Learning Large Incomplete Matrices. \emph{Journal of Machine 
#' Learning Research}, \bold{11}(80), 2287--2322.
#' 
#' @seealso \code{\link{soft_impute_tune}()}, \code{\link{fraction_grid}()}
#' 
#' @examples
#' # toy example derived from MovieLens 100K dataset
#' data("MovieLensToy")
#' # Soft-Impute with discretization step
#' fit <- soft_impute(MovieLensToy, discretize = TRUE)
#' # extract discretized completed matrix with fifth value 
#' # of regularization parameter
#' X_hat <- get_completed(fit, which = 5, discretized = TRUE)
#' head(X_hat)
#' 
#' @keywords multivariate
#' 
#' @importFrom softImpute biScale complete softImpute
#' @export

soft_impute <- function(X, lambda = fraction_grid(reverse = TRUE), 
                        relative = TRUE, type = c("svd", "als"), 
                        rank.max = NULL, thresh = 1e-05, 
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
  if (isTRUE(discretize)) {
    # if not supplied, obtain values of rating scale (categories)
    if (is.null(values)) values <- unique(X[!is.na(X)])
    # obtain minimum and maximum value for discretization
    min_value <- min(values)
    max_value <- max(values)
    # discretize the imputed matrix
    # FIXME: this only works for integer sequence from min_value to max_value
    if (nb_lambda == 1L) {
      out$X_discretized <- pmin(pmax(round(X_imputed), min_value), max_value)
    } else {
      out$X_discretized <- lapply(X_imputed, function(current_X_imputed) {
        pmin(pmax(round(current_X_imputed), min_value), max_value)
      })
    }
  }
  
  # add class and return object
  class(out) <- "soft_impute"
  out
}
