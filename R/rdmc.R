# ************************************
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ************************************


#' @useDynLib rdmc, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @export

rdmc <- function(X, values = NULL, lambda = fraction_grid(), relative = TRUE, 
                 rank_max = NULL, type = "svd", svd_tol = 1e-05, 
                 loss = c("pseudo_huber", "absolute", "bounded"),
                 loss_const = NULL, delta = 1.05, mu = 0.1, 
                 conv_tol = 1e-02, max_iter = 10L, 
                 # starting values
                 L = NULL, Theta = NULL) {
  
  # initializations
  X <- as.matrix(X)
  d <- dim(X)
  
  # construct indicator matrix of missing values
  is_NA <- is.na(X)
  
  # check values of rating scale
  if (is.null(values)) values <- unique(X[!is_NA])
  values <- sort(values)  # ensure values of rating scale are sorted
  # check values of tuning parameter
  lambda <- sort(unique(lambda))  # ensure values of tuning parameter are sorted
  relative <- isTRUE(relative)
  # check maximum rank
  if (is.null(rank_max)) rank_max <- min(dim(X))
  # check loss function
  loss <- match.arg(loss)
  if (is.null(loss_const)) {
    # set default constant for loss function (if applicable)
    loss_const <- switch(loss, pseudo_huber = 1, absolute = NA_real_, 
                         bounded = (max(values) - min(values)) / 2)
  }
  
  # center data matrix with columnwise median of observed data
  medians <- apply(X, 2L, median, na.rm = TRUE)
  X <- sweep(X, 2, medians, FUN = "-")
  # update discrete candidate values for each column
  values <- sapply(medians, function(m) values - m)
  
  # if relative grid of tuning parameter values is requested, compute the 
  # largest singular value of the median-imputed matrix
  if (relative) {
    X_zeros <- X
    X_zeros[is_NA] <- 0
    d_max <- svd(X_zeros)$d[1L]  # largest singular value
  } else d_max <- 1
  
  # initialize starting values if not supplied
  # (only use user-supplied starting values if both are supplied)
  if (is.null(L) || is.null(Theta)) {
    # since X is median-centered, we can initialize missing values in L with 0
    if (relative) L <- X_zeros
    else {
      L <- X
      L[is_NA] <- 0
    }
    # initialize Theta with 0
    Theta <- matrix(0, nrow = d[1], ncol = d[2])
  }

  # convert indicator matrix of missing values to matrix of row and column 
  # indices of missing and observed elements
  idx_NA <- which(is_NA, arr.ind = TRUE, useNames = FALSE)
  idx_observed <- which(!is_NA, arr.ind = TRUE, useNames = FALSE)
  
  # call C++ function
  out <- rdmc_cpp(X, idx_NA = idx_NA - 1L, idx_observed = idx_observed - 1L, 
                  values = values, lambda = lambda, d_max = d_max, 
                  rank_max = rank_max, type = type, svd_tol = svd_tol, 
                  loss = loss, loss_const = loss_const, delta = delta, 
                  mu = mu, conv_tol = conv_tol, max_iter = max_iter,  
                  L = L, Theta = Theta)
  
  # obtain completed matrix on original rating scale
  if (length(lambda) == 1L) {
    X[is_NA] <- out$L[is_NA]
    out$X <- sweep(X, 2L, medians, FUN = "+")
  } else {
    out$X <- lapply(out$L, function(L) {
      X[is_NA] <- L[is_NA]
      sweep(X, 2L, medians, FUN = "+")
    })
  }
  
  # add class and return object
  class(out) <- "rdmc"
  out
  
}
