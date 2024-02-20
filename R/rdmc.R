# **************************************
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# **************************************

#' @useDynLib rdmc, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @export

rdmc <- function(X, values, lambda, type = "svd", loss = "pseudo_huber", 
                 svd_tol = 1e-05, delta = 1.05, mu = 0.1, conv_tol = 1e-02, 
                 max_iter = 10) {
  
  # check arguments
  values <- sort(values)       # make sure values of rating scale are sorted
  nb_values <- length(values)  # number of values in rating scale (categories)
  type <- match.arg(type)
  
  # construct indicator matrix of missing values (as integers)
  is_NA <- is.na(X)
  storage.mode(is_NA) <- "integer"

  # center data matrix with midpoint of rating scale
  midpoint <- mean(values[c(1, nb_values)])
  X <- X - midpoint
  values <- values - midpoint
  
  # call C++ function
  out <- rdmc_cpp(X, is_NA = is_NA, values = values, lambda = lambda, 
                  type = type, svd_tol = svd_tol, delta = delta, mu = mu, 
                  conv_tol = conv_tol, max_iter = max_iter)
  
  # obtain completed matrix on original rating scale
  storage.mode(is_NA) <- "logical"
  if (length(lambda) == 1L) {
    X[is_NA] <- out$L[is_NA]
    out$X <- X + midpoint
  } else {
    # restructure output from C++
    out <- list(
      lambda = sapply(out, "[[", "lambda"),
      L = lapply(out, "[[", "L"),
      Z = lapply(out, "[[", "Z"),
      Theta = lapply(out, "[[", "Theta"),
      objective = sapply(out, "[[", "objective"),
      converged = sapply(out, "[[", "converged"),
      nb_iter = sapply(out, "[[", "nb_iter")
    )
    # add completed matrix on original rating scale
    out$X <- lapply(out$L, function(L) {
      X[is_NA] <- L[is_NA]
      X + midpoint
    })
  }
  
  # add class and return object
  class(out) <- "rdmc"
  out
}
