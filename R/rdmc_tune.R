# ************************************
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ************************************


#' Robust discrete matrix completion with hyperparameter tuning
#' 
#' Perform robust discrete matrix completion with a low-rank constraint on a 
#' latent continuous matrix, implemented via an ADMM algorithm.  The 
#' regularization parameter is thereby selected via repeated holdout validation 
#' or cross-validation.
#' 
#' @inherit rdmc details
#' 
#' @inheritParams rdmc
#' @param splits  an object inheriting from class \code{"split_control"}, as 
#' generated by \code{\link{holdout_control}()} for repeated holdout validation 
#' or \code{\link{cv_folds_control}()} for \eqn{K}-fold cross-validation, or a 
#' list of index vectors giving different validation sets of observed cells as 
#' generated by \code{\link{create_splits}()}.  Cells in the validation set
#' will be set to \code{NA} for fitting the algorithm with the training set of 
#' observed cells.
#' @param \dots  additional arguments to be passed down to \code{\link{rdmc}()}.
#' 
#' @return 
#' An object of class \code{"rdmc_tuned"} with the following components: 
#' \item{lambda}{a numeric vector containing the values of the regularization 
#' parameter.}
#' \item{tuning_loss}{a numeric vector containing the (average) values of the 
#' loss function on the validation set(s) for each value of the regularization 
#' parameter.}
#' \item{lambda_opt}{numeric; the optimal value of the regularization 
#' parameter.}
#' \item{fit}{an object of class \code{"\link{rdmc}"} containing the results 
#' from the algorithm with the optimal regularization parameter on the full 
#' (observed) data matrix.}
#' 
#' The class structure is still experimental and may change in the future. 
#' The following accessor functions are available:
#' \itemize{
#'   \item \code{\link{get_completed}()} to extract the completed (i.e., 
#'   imputed) data matrix with the optimal value of the regularization 
#'   parameter,
#'   \item \code{\link{get_lambda}()} to extract the optimal value of the 
#'   regularization parameter,
#'   \item \code{\link{get_nb_iter}()} to extract the number of iterations with 
#'   the optimal value of the regularization parameter.
#' }
#' 
#' @author Andreas Alfons
#' 
#' @inherit rdmc references
#' 
#' @seealso 
#' \code{\link{rdmc}()}, \code{\link{fraction_grid}()}, 
#' 
#' \code{\link{holdout_control}()}, \code{\link{cv_folds_control}()}, 
#' \code{\link{create_splits}()}
#' 
#' @examples
#' # toy example derived from MovieLens 100K dataset
#' data("MovieLensToy")
#' # robust discrete matrix completion with hyperparameter tuning
#' set.seed(20250723)
#' fit <- rdmc_tune(MovieLensToy, 
#'                  lambda = fraction_grid(nb_lambda = 6),
#'                  splits = holdout_control(R = 5))
#' # extract completed matrix with optimal regularization parameter
#' X_hat <- get_completed(fit)
#' head(X_hat)
#' # extract optimal value of regularization parameter
#' get_lambda(fit)
#' # extract number of iterations with optimal regularization parameter
#' get_nb_iter(fit)
#' 
#' @keywords multivariate
#' 
#' @export

rdmc_tune <- function(X, values = NULL, lambda = fraction_grid(), 
                      relative = TRUE, splits = holdout_control(), 
                      loss = c("pseudo_huber", "absolute", "truncated"),
                      loss_const = NULL, ...) {
  
  # initializations
  X <- as.matrix(X)
  
  # construct index vector of observed values
  observed <- which(!is.na(X))
  
  # check values of rating scale
  if (is.null(values)) values <- unique(X[observed])
  values <- sort(values)  # ensure values of rating scale are sorted
  # check values of tuning parameter
  lambda <- sort(unique(lambda))  # ensure values of tuning parameter are sorted
  if (length(lambda) == 1L) {
    stop("only one value of 'lambda'; use function rdmc() instead")
  }
  relative <- isTRUE(relative)
  # check loss function
  loss <- match.arg(loss)
  if (is.null(loss_const)) {
    # set default constant for loss function (if applicable)
    loss_const <- switch(
      loss, 
      pseudo_huber = (max(values) - min(values)) / (length(values) - 1), 
      absolute = NA_real_, 
      truncated = (max(values) - min(values)) / 2
    )
  }
  
  # create splits for tuning parameter validation
  if (inherits(splits, "split_control")) {
    splits <- create_splits(observed, control = splits)
  }
  
  # fit robust discrete matrix completion to the different training data sets
  fit_train <- lapply(splits, function(indices, ...) {
    # create training data where elements from test set are set to NA
    X_train <- X
    X_train[indices] <- NA_real_
    # apply robust discrete matrix completion to training data
    fit_train <- rdmc(X_train, values = values, lambda = lambda, 
                      relative = relative, loss = loss, 
                      loss_const = loss_const, ...)
  }, ...)
  
  # extract predictions for the elements in the different test sets and compute 
  # prediction loss
  tuning_loss <- mapply(function(indices, fit) {
    # extract elements from test set as a vector
    X_test <- X[indices]
    # compute prediction loss
    sapply(fit$X, function(X_hat) {
      if (loss == "pseudo_huber") {
        pseudo_huber(X_test - X_hat[indices], const = loss_const)
      } else if (loss == "absolute") {
        abs(X_test - X_hat[indices])
      } else {
        truncated(X_test - X_hat[indices], const = loss_const)
      }
    })
    
  }, indices = splits, fit = fit_train, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  # combine results for tuning parameter validation into one matrix
  # (each column holds the prediction losses for the corresponding lambda)
  tuning_loss <- do.call(rbind, tuning_loss)
  # compute column means
  tuning_loss <- colMeans(tuning_loss)
  
  # select the optimal lambda: reverse vectors so that in the unlikely case of 
  # ties, we select the lambda with stronger penalization
  which_opt <- rev(seq_along(lambda))[which.min(rev(tuning_loss))]
  lambda_opt <- lambda[which_opt]
  
  # Note: It's possible that on different training sets, we get different 
  # medians in some variables so that the discrete constraint is different 
  # for those variables. But for the starting values, it shouldn't matter 
  # (L doesn't even have to satisfy any discrete constraint), as we simply 
  # use L and Theta for the first soft-thresholding SVD step to update Z, 
  # and subsequently L is updated to satisfy the discrete constrains using 
  # the medians on the full data set. Hence we could just take the mean or 
  # median over the training sets for L rather than the mode (which currently 
  # uses random sampling in case of multiple modes). Specifically, we apply 
  # the soft-thresholding SVD step to L + Theta/mu, so it may be smoother to 
  # use the linearity of the mean to compute both starting values for L and 
  # Theta, rather than destroying the linear relationship by using the mode 
  # for L and the mean for Theta.
  
  # compute starting values for L and Theta by aggregating solutions obtained 
  # from the different training data
  L <- sapply(fit_train, function(fit) fit$L[[which_opt]], 
              simplify = "array")
  L <- apply(L, 1:2, median)
  Theta <- sapply(fit_train, function(fit) fit$Theta[[which_opt]], 
                  simplify = "array")
  Theta <- apply(Theta, 1:2, mean)
  
  # apply robust discrete matrix completion with optimal tuning parameter
  fit_opt <- rdmc(X, values = values, lambda = lambda_opt, relative = relative, 
                  loss = loss, loss_const = loss_const, ..., L = L, 
                  Theta = Theta)
  
  # construct list of relevant output
  out <- list(lambda = lambda, tuning_loss = tuning_loss, 
              lambda_opt = lambda_opt, fit = fit_opt)
  class(out) <- "rdmc_tuned"
  
  # return output
  out
  
}


# ## function to compute mode for starting value of L in fit with optimal lambda
# randomized_mode <- function(x) {
#   # compute contingency table
#   tab <- table(x)
#   # determine the maximum frequency
#   max <- max(tab)
#   # keep all values that occur with the maximum frequency
#   out <- as.numeric(names(tab)[tab == max])
#   # if there are multiple modes, randomly select one to be returned
#   if (length(out) > 1) sample(out, 1)
#   else out
# }
