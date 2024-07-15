# ************************************
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ************************************

#' @useDynLib rdmc, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @export

rdmc <- function(X, values = NULL, lambda, type = "svd", 
                 loss = c("pseudo_huber", "absolute", "bounded"), 
                 svd_tol = 1e-05, loss_tuning = NULL, delta = 1.05, 
                 mu = 0.1, conv_tol = 1e-02, max_iter = 10, 
                 # starting values: only used for fitting the algorithm for 
                 # the optimal lambda after tuning, but they are currently 
                 # ignored if a vector of lambda is supplied
                 L = NULL, Z = NULL, Theta = NULL) {
  
  # initializations
  X <- as.matrix(X)
  
  # construct indicator matrix of missing values
  is_NA <- is.na(X)
  
  # check arguments
  if (is.null(values)) values <- unique(X[!is_NA])
  values <- sort(values)         # ensure values of rating scale are sorted
  nb_values <- length(values)    # number of values in rating scale (categories)
  lambda <- sort(unique(lambda)) # ensure values of tuning parameter are sorted
  type <- match.arg(type)
  loss <- match.arg(loss)
  
  # check bound in case bounded loss function
  if (is.null(loss_tuning)) {
    loss_tuning <- switch(loss, 
                          pseudo_huber = 1,
                          absolute = NA_real_,
                          bounded = (max(values) - min(values)) / 2)
  }
  
  # # center data matrix with midpoint of rating scale
  # midpoint <- mean(values[c(1, nb_values)])
  # X <- X - midpoint
  # values <- values - midpoint
  # center data with columnwise median of observed data
  # TODO: there is potential for optimization, as this can be computed in 
  # rdmc_tune() and passed down. In addition, since we compute the medians 
  # already in R, they can be passed down to C++ to initialize the starting 
  # values for the first value of lambda.
  medians <- apply(X, 2, median, na.rm = TRUE)
  X <- sweep(X, 2, medians, FUN = "-")
  # update discrete candidate values for each column
  value_mat <- sapply(medians, function(m) values - m)
  
  # call C++ function (requires indicator matrix as integers)
  storage.mode(is_NA) <- "integer"
  if (length(lambda) == 1L && !is.null(L) & !is.null(Z) && !is.null(Theta)) {
    out <- rdmc_opt_cpp(X, is_NA = is_NA, values = value_mat, lambda = lambda, 
                        type = type, loss = loss, svd_tol = svd_tol, 
                        loss_tuning = loss_tuning, delta = delta, mu = mu, 
                        conv_tol = conv_tol, max_iter = max_iter, L = L, 
                        Z = Z, Theta = Theta)
  } else {
    out <- rdmc_cpp(X, is_NA = is_NA, values = value_mat, lambda = lambda, 
                    type = type, loss = loss, svd_tol = svd_tol, 
                    loss_tuning = loss_tuning, delta = delta, mu = mu, 
                    conv_tol = conv_tol, max_iter = max_iter)
  }
  
  # obtain completed matrix on original rating scale
  storage.mode(is_NA) <- "logical"
  if (length(lambda) == 1L) {
    X[is_NA] <- out$L[is_NA]
    # out$X <- X + midpoint
    out$X <- sweep(X, 2, medians, FUN = "+")
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
      # X + midpoint
      sweep(X, 2, medians, FUN = "+")
    })
  }
  
  # add class and return object
  class(out) <- "rdmc"
  out
}
