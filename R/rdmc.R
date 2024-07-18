# ************************************
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ************************************

#' @useDynLib rdmc, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @export

rdmc <- function(X, values = NULL, lambda, rank_max = NULL, 
                 type = "svd", svd_tol = 1e-05, 
                 loss = c("bounded", "absolute", "pseudo_huber"),
                 loss_const = NULL, delta = 1.05, mu = 0.1, 
                 conv_tol = 1e-02, max_iter = 10, 
                 # starting values
                 L = NULL, Theta = NULL) {
  
  # initializations
  X <- as.matrix(X)
  d <- dim(X)
  
  # construct indicator matrix of missing values
  is_NA <- is.na(X)
  
  # check arguments
  if (is.null(values)) values <- unique(X[!is_NA])
  values <- sort(values)          # ensure values of rating scale are sorted
  lambda <- sort(unique(lambda))  # ensure values of tuning parameter are sorted
  if (is.null(rank_max)) rank_max <- min(dim(X))
  loss <- match.arg(loss)
  if (is.null(loss_const)) {
    # set default constant for loss function (if applicable)
    loss_const <- switch(loss, bounded = (max(values) - min(values)) / 2, 
                         absolute = NA_real_, pseudo_huber = 1)
  }
  
  # center data matrix with columnwise median of observed data
  medians <- apply(X, 2L, median, na.rm = TRUE)
  X <- sweep(X, 2, medians, FUN = "-")
  # update discrete candidate values for each column
  values <- sapply(medians, function(m) values - m)
  
  # initialize starting values if not supplied
  # (only use user-supplied starting values if both are supplied)
  if (is.null(L) || is.null(Theta)) {
    # since X is median-centered, we can initialize missing values in L with 0
    L <- X
    L[is_NA] <- 0
    # initialize Theta with 0
    Theta <- matrix(0, nrow = d[1], ncol = d[2])
  }

  # convert indicator matrix of missing values to matrix of row and column 
  # indices of missing and observed elements
  idx_NA <- which(is_NA, arr.ind = TRUE, useNames = FALSE)
  idx_observed <- which(!is_NA, arr.ind = TRUE, useNames = FALSE)
  
  # call C++ function
  out <- rdmc_cpp(X, idx_NA = idx_NA - 1L, idx_observed = idx_observed - 1L, 
                  values = values, lambda = lambda, rank_max = rank_max, 
                  type = type, svd_tol = svd_tol, loss = loss, 
                  loss_const = loss_const, delta = delta, mu = mu, 
                  conv_tol = conv_tol, max_iter = max_iter,  L = L, 
                  Theta = Theta)

  # obtain completed matrix on original rating scale
  if (length(lambda) == 1L) {
    X[is_NA] <- out$L[is_NA]
    out$X <- sweep(X, 2L, medians, FUN = "+")
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
      sweep(X, 2L, medians, FUN = "+")
    })
  }
  
  # add class and return object
  class(out) <- "rdmc"
  out
}
