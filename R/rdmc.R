# **************************************
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# **************************************

#' @useDynLib rdmc, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @export

rdmc <- function(X, nb_cat, lambda, type = "svd", loss = "pseudo_huber", 
                 svd_tol = 1e-05, delta = 1.05, mu = 0.1, conv_tol = 1e-02, 
                 max_iter = 10) {
  
  # check arguments
  type <- match.arg(type)
  
  # construct indicator matrix of missing values (as integers)
  is_NA <- is.na(X)
  storage.mode(is_NA) <- "integer"

  # call C++ function
  rdmc_cpp(X, is_NA = is_NA, nb_cat = nb_cat, lambda = lambda, type = type, 
           svd_tol = svd_tol, delta = delta, mu = mu, conv_tol = conv_tol, 
           max_iter = max_iter)
}
