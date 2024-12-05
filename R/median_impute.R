# ************************************
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ************************************

#' Median imputation
#' 
#' Perform median imputation.  In case of discrete rating-scale data, a 
#' discretization step can be carried out afterwards to make sure that the 
#' imputed values are mapped to the rating scale of the observed values (as 
#' the median of a given column may lie in between two answer categories in 
#' case of an even number of observed values).  This is done by randomly 
#' sampling from the largest answer category smaller than the median and the 
#' smallest answer category larger than the median (for each missing cell).
#' 
#' @param X  a matrix or data frame with missing values.
#' @param discretize  a logical indicating whether to include a discretization 
#' step after median imputation (defaults to \code{TRUE}).  In case of 
#' discrete rating-scale data, this can be used to ensure that the imputed 
#' values are mapped to the discrete rating scale of the observed values.
#' @param values  an optional numeric vector giving the possible values of 
#' discrete ratings.  This is ignored if \code{discretize} is \code{FALSE}.  
#' Currently, the possible values are assumed to be the same for all columns.  
#' If \code{NULL}, the unique values of the observed parts of \code{X} are 
#' used.
#' 
#' @return 
#' An object of class \code{"median_impute"}.
#' 
#' The class structure is still experimental and may change.  Use the accessor 
#' function \code{\link{get_X}()} to extract the imputed data matrix.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{mode_impute}()}
#' 
#' @keywords multivariate
#' 
#' @importFrom stats median
#' @export
median_impute <- function(X, discretize = TRUE, values = NULL) {
  
  # initializations
  X <- as.matrix(X)
  p <- ncol(X)
  seq_p <- seq_len(p)
  discretize <- isTRUE(discretize)
  
  # construct indicator matrix of missing values
  is_NA <- is.na(X)
  
  # compute columnwise medians
  medians <- apply(X, 2, median, na.rm = TRUE)
  
  # impute with columnwise medians
  X_imputed <- lapply(seq_p, function(j) {
    X_j <- X[, j]
    is_NA_j <- is_NA[, j]
    X_j[is_NA_j] <- medians[j]
    X_j
  })
  X_imputed <- do.call(cbind, X_imputed)
  colnames(X_imputed) <- colnames(X)
  
  # construct object to be returned
  out <- list(medians = medians, X = X_imputed)
  
  # if requested, discretize the imputed matrix in case the median falls in 
  # between answer categories
  if (discretize) {
    # check values of rating scale
    if (is.null(values)) values <- unique(X[!is_NA])
    values <- sort(values)  # ensure values of rating scale are sorted
    # find which columns need to be discretized
    discretize <- !(medians %in% values)
    # discretize the imputed matrix
    X_discretized <- lapply(seq_p, function(j) {
      X_j <- X_imputed[, j]
      if (discretize[j]) {
        # obtain number of missing values
        is_NA_j <- is_NA[, j]
        n_NA_j <- sum(is_NA_j)
        # find largest answer category smaller than median and smallest 
        # answer category larger than median
        which_lower <- rev(which(values < medians[j]))[1L]
        which_upper <- which(values > medians[j])[1L]
        candidates <- values[c(which_lower, which_upper)]
        # replace missing values with random draw from those answer categories
        X_j[is_NA_j] <- sample(candidates, n_NA_j, replace = TRUE)
      }
      X_j
    })
    X_discretized <- do.call(cbind, X_discretized)
    colnames(X_discretized) <- colnames(X)
    # add discretized matrix to object to be returned
    out$X_discretized <- X_discretized
  }
  
  # add class and return object
  class(out) <- "median_impute"
  out
  
}
