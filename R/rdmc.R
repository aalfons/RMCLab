# ************************************
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ************************************


#' Robust discrete matrix completion
#' 
#' Perform robust discrete matrix completion with a low-rank constraint on a 
#' latent continuous matrix, implemented via an ADMM algorithm.
#' 
#' For the loss part of the objective function, the pseudo-Huber loss 
#' (\code{loss = "pseudo_huber"}) is given by
#' \deqn{\rho(x) = \code{loss\_const}^2 (\sqrt{1 + (x/\code{loss\_const})^2} - 1).}{rho(x) = \code{loss_const}^2 (sqrt(1 + (x/\code{loss_const})^2) - 1).}
#' The absolute loss
#' (\code{loss = "absolute"}) is given by
#' \deqn{\rho(x) = |x|,}{rho(x) = |x|,}
#' and the bounded absolute loss (\code{loss = "bounded"}) is defined as
#' \deqn{\rho(x) = \min (|x|, \code{loss\_const}).}{rho(x) = min(|x|, \code{loss_const}).}
#' 
#' @param X  a matrix or data frame of discrete ratings with missing values.
#' @param values  an optional numeric vector giving the possible values of the 
#' ratings.  Currently, these are assumed to be the same for all columns.  If 
#' \code{NULL}, the unique values of the observed parts of \code{X} are used.
#' @param lambda  a numeric vector giving values of the regularization 
#' parameter.  See \code{\link{fraction_grid}()} for the default values.
#' @param relative  a logical indicating whether the values of the 
#' regularization parameter should be considered relative to a certain 
#' reference value computed from the data at hand.  If \code{TRUE} (the 
#' default), the values of \code{lambda} are multiplied with the largest 
#' singular value of the median-centered data matrix with missing values 
#' replaced by zeros.
#' @param type  a character string specifying the type of algorithm for the 
#' low-rank latent continuous matrix.  Currently only \code{"svd"} is 
#' implemented for a soft-thresholded SVD step.
#' @param loss  a character string specifying the robust loss function for the 
#' loss part of the objective function.  Possible values are 
#' \code{"pseudo_huber"} (the default) for the pseudo-Huber loss, 
#' \code{"absolute"} for the absolute loss, and \code{"bounded"} for the 
#' bounded absolute loss.  See \sQuote{Details} for more information.
#' @param loss_const  tuning constant for the loss function.  For the 
#' pseudo-Huber loss, the default value is the average step size between the 
#' rating categories in \code{values}.  For the bounded absolute loss, 
#' the default is half the range of the rating categories in \code{values}.  
#' This is ignored for the absolute loss, which does not have a tuning 
#' parameter.  See \sQuote{Details} for more information.
#' @param svd_tol  numeric tolerance for the soft-thresholded SVD step.  Only 
#' singular values larger than \code{svd_tol} are kept to construct the 
#' low-rank latent continuous matrix.
#' @param rank_max  a positive integer giving a rank constraint in the 
#' soft-thresholded SVD step for the latent continuous matrix. The default is 
#' to use the minimum of the number of rows and columns.
#' @param mu  numeric; penalty parameter for the discrepancy between the 
#' discrete rating matrix and the latent low-rank continuous matrix.  It is 
#' not recommended to change the default value of 0.1.
#' @param delta  numeric; update factor for penalty parameter \code{mu} applied 
#' after each iteration to increase the strength of the penalty.  It is not 
#' recommended to change the default value of 1.05.
#' @param conv_tol  numeric; convergence tolerance for the relative change in 
#' the objective function.
#' @param max_iter  a positive integer specifying the maximum number of 
#' iterations.  In practice, large gains can often be had in the first few 
#' iterations, with subsequent iterations yielding relatively small gains until 
#' convergence.  Hence the default is to perform at most 10 iterations.
#' @param L,Theta  starting values for the algorithm.  These are not expected 
#' to be set by the user.  Instead, it is recommended to call this function 
#' with a grid of values for the regularization parameter \code{lambda} so that 
#' the implementation automatically takes advantage of warm starts.
#' 
#' @return 
#' An object of class \code{"rdmc"}.  The class structure is still experimental 
#' and may change.  
#' 
#' @author Andreas Alfons and Aurore Archimbaud
#' 
#' @references
#' Archimbaud, A., Alfons, A., and Wilms, I. (2024) Robust Matrix Completion 
#' for Discrete Rating-Scale Data. arXiv:2412.20802. 
#' \doi{10.48550/arXiv.2412.20802}.
#' 
#' @seealso \code{\link{rdmc_tune}()}, \code{\link{fraction_grid}()}
#' 
#' @keywords multivariate
#' 
#' @useDynLib rdmc, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @export

rdmc <- function(X, values = NULL, lambda = fraction_grid(), relative = TRUE, 
                 loss = c("pseudo_huber", "absolute", "bounded"),
                 loss_const = NULL, type = "svd", svd_tol = 1e-04, 
                 rank_max = NULL, mu = 0.1, delta = 1.05, 
                 conv_tol = 1e-04, max_iter = 100L, 
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
    loss_const <- switch(
      loss, 
      pseudo_huber = (max(values) - min(values)) / (length(values) - 1), 
      absolute = NA_real_, 
      bounded = (max(values) - min(values)) / 2
    )
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
