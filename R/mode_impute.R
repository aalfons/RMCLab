# ************************************
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ************************************


#' Mode imputation
#' 
#' Perform mode imputation for discrete data.  In case of multiple modes in a 
#' given column, one of them is selected at random for each missing cell.
#' 
#' @param X  a matrix or data frame of discrete data with missing values.
#' @param values  an optional numeric vector giving the possible values.  
#' Currently, the possible values are assumed to be the same for all columns.  
#' If \code{NULL}, the unique values of the observed parts of \code{X} are 
#' used.
#' 
#' @return 
#' An object of class \code{"mode_impute"}.
#' 
#' The class structure is still experimental and may change.  Use the accessor 
#' function \code{\link{get_X}()} to extract the imputed data matrix.
#' 
#' @note  
#' The mode is computed as the most frequent value, hence this function is only 
#' suitable for discrete data.  It does not estimate the mode of a continuous 
#' density.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{median_impute}()}
#' 
#' @keywords multivariate
#' 
#' @export

mode_impute <- function(X, values = NULL) {
  
  # initializations
  X <- as.matrix(X)
  p <- ncol(X)
  seq_p <- seq_len(p)
  
  # construct indicator matrix of missing values
  is_NA <- is.na(X)
  
  # check values of rating scale
  if (is.null(values)) values <- unique(X[!is_NA])
  values <- sort(values)  # ensure values of rating scale are sorted
  
  # compute columnwise medians
  modes <- apply(X, 2, .mode, values = values, simplify = FALSE)
  
  # impute with columnwise medians
  X_imputed <- lapply(seq_p, function(j) {
    X_j <- X[, j]
    is_NA_j <- is_NA[, j]
    modes_j <- modes[[j]]
    if (length(modes_j) == 1L) {
      # single mode
      X_j[is_NA_j] <- modes_j
    } else {
      # replace missing values with random draw from multiple modes
      n_NA_j <- sum(is_NA_j)
      X_j[is_NA_j] <- sample(modes_j, n_NA_j, replace = TRUE)
    }
    X_j
  })
  X_imputed <- do.call(cbind, X_imputed)
  colnames(X_imputed) <- colnames(X)
  
  # construct object to be returned
  out <- list(modes = modes, X = X_imputed)
  class(out) <- "mode_impute"
  out
  
}


## internal function to compute mode
.mode <- function(x, values) {
  # compute contingency table
  tab <- tabulate(factor(x, levels = values, exclude = NA), 
                  nbins = length(values))
  # determine the maximum frequency
  max <- max(tab)
  # keep all values that occur with the maximum frequency
  values[tab == max]
}
