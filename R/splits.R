# ************************************
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ************************************


#' @name validation_control
#'
#' @title Control objects for hyperparameter validation
#' 
#' @description Construct control objects for repeated holdout validation or 
#' \eqn{K}-fold cross-validation.
#' 
#' @return An object inheriting from class \code{"split_control"} containing 
#' the relevant information for splitting the the observed cells of a data 
#' matrix into training and validation sets for hyperparameter tuning.
#' 
#' The subclass \code{"holdout_control"} returned by \code{holdout_control()} 
#' is a list with components \code{pct} and \code{R} containing the 
#' corresponding argument values after validity checks.
#' 
#' The subclass \code{"cv_folds_control"} returned by \code{cv_folds_control()} 
#' is a list with a single component \code{K} containing the corresponding 
#' argument value after validity checks.
#' 
#' @seealso 
#' \code{\link{create_splits}()}, 
#' 
#' \code{\link{rdmc_tune}()}, \code{\link{soft_impute_tune}()}
#' 
#' @examples
#' # toy example derived from MovieLens 100K dataset
#' data("MovieLensToy")
#' # robust discrete matrix completion with hyperparameter tuning
#' set.seed(20250723)
#' fit <- rdmc_tune(MovieLensToy, 
#'                  lambda = fraction_grid(nb_lambda = 6),
#'                  splits = holdout_control(R = 5))
#' # extract optimal value of regularization parameter
#' get_lambda(fit)
#' 
#' @keywords utilities

NULL


## control object for holdout validation
#' @rdname validation_control
#' 
#' @param pct  numeric in the interval (0, 1); the percentage of observed cells 
#' in the data matrix to be randomly selected into the validation set (defaults 
#' to 0.1).
#' @param R  an integer giving the number of random splits into training and 
#' validation sets (defaults to 10).
#' 
#' @export
holdout_control <- function(pct = 0.1, R = 10L) {
  if (pct <= 0 || pct >= 1) {
    stop("the percentage of holdout elements must be in the interval (0, 1)")
  }
  if (R < 1L) stop("holdout validation requires at least 1 replication")
  control <- list(pct = pct, R = R)
  class(control) <- c("holdout_control", "split_control")
  control
}


## control object for cross-validation
#' @rdname validation_control
#' 
#' @param K  an integer giving the number of cross-validation folds (defaults 
#' to 5).
#' 
#' @export
cv_folds_control <- function(K = 5L) {
  if (K < 2L) stop("cross-validation requires at least 2 folds")
  control <- list(K = K)
  class(control) <- c("cv_folds_control", "split_control")
  control
}


## generic function for validation splits
#' Create splits of observed data cells for hyperparameter tuning
#' 
#' Split the observed cells of a data matrix into training and validation sets 
#' for hyperparameter tuning.  Methods are available for repeated holdout 
#' validation and \eqn{K}-fold cross-validation.
#' 
#' Functions \code{holdout()} and \code{cv_folds()} are wrapper 
#' functions that first call \code{\link{holdout_control}()} and 
#' \code{\link{cv_folds_control}()}, respectively, before calling 
#' \code{create_splits()}.
#' 
#' @param indices  an integer vector giving the indices of observed cells in a 
#' data matrix.
#' @param control  a control object inheriting from class 
#' \code{"split_control"} as generated by \code{\link{holdout_control}()} 
#' for repeated holdout validation or \code{\link{cv_folds_control}()} for 
#' \eqn{K}-fold cross-validation.
#' @inheritParams validation_control
#' 
#' @return  A list of index vectors giving the validation sets of the 
#' respective replication or cross-validation fold.
#' 
#' @author Andreas Alfons
#' 
#' @seealso 
#' \code{\link{holdout_control}()}, \code{\link{cv_folds_control}()}, 
#' 
#' \code{\link{rdmc_tune}()}, \code{\link{soft_impute_tune}()}
#' 
#' @examples
#' # toy example derived from MovieLens 100K dataset
#' data("MovieLensToy")
#' # set up validation sets so that methods use same data splits
#' set.seed(20250723)
#' observed <- which(!is.na(MovieLensToy))
#' holdout_splits <- holdout(observed, R = 5)
#' # robust discrete matrix completion with hyperparameter tuning
#' fit_RDMC <- rdmc_tune(
#'   MovieLensToy, 
#'   lambda = fraction_grid(nb_lambda = 6),
#'   splits = holdout_splits
#' )
#' # Soft-Impute with discretization step and hyperparameter tuning
#' fit_SI <- soft_impute_tune(
#'   MovieLensToy, 
#'   lambda = fraction_grid(nb_lambda = 6, reverse = TRUE),
#'   splits = holdout_splits
#' )
#' # extract optimal values of regularization parameter
#' get_lambda(fit_RDMC)
#' get_lambda(fit_SI)
#' 
#' @keywords utilities
#' 
#' @export

create_splits <- function(indices, control) UseMethod("create_splits", control)


## method to draw indices for holdout validation
#' @export
create_splits.holdout_control <- function(indices, control) {
  # initializations
  n <- length(indices)
  pct <- control$pct
  R <- control$R
  # check if enough indices have been supplied
  min <- ceiling(0.5000000001/pct) # avoid issues with floating point arithmetic
  if (n < min) {
    stop(sprintf("'indices' should contain at least %d elements", min))
  }
  # draw indices of holdout set
  m <- round(pct * n)
  replicate(R, sample(indices, size = m), simplify = FALSE)
}


## method to draw folds of indices for K-fold cross-validation
#' @export
create_splits.cv_folds_control <- function(indices, control) {
  # initialization
  n <- length(indices)
  K <- control$K
  # check if enough indices have been supplied
  if (n < K) stop(sprintf("'indices' should contain at least %d elements", K))
  # permute observations
  indices <- sample(indices)
  # assign a block to each observation
  blocks <- rep(seq_len(K), length.out = n)
  # split the permuted observations according to the block they belong to
  folds <- split(indices, blocks)
  names(folds) <- NULL
  folds
}
  

## wrapper function to draw indices for holdout validation
#' @rdname create_splits
#' @export
holdout <- function(indices, pct = 0.1, R = 10L) {
  control <- holdout_control(pct = pct, R = R)
  create_splits(indices, control)
}


## wrapper function to draw folds of indices for K-fold cross-validation
#' @rdname create_splits
#' @export
cv_folds <- function(indices, K = 5L) {
  control <- cv_folds_control(K = K)
  create_splits(indices, control)
}


# ## function to draw indices for holdout validation
# holdout <- function(indices, pct = 0.1, R = 10) {
#   n <- length(indices)
#   min <- ceiling(0.5000000001/pct) # avoid issues with floating point arithmetic
#   if (n < min) {
#     stop(sprintf("'indices' should contain at least %d elements", min))
#   }
#   m <- round(pct * n)
#   replicate(R, sample(indices, size = m), simplify = FALSE)
# }
# 
# ## function to draw folds of indices for K-fold cross-validation
# cv_folds <- function(indices, K = 5L) {
#   n <- length(indices)
#   if (K < 2L) stop("cross-validation requires at least 2 folds")
#   if (n < K) stop(sprintf("'indices' should contain at least %d elements", K))
#   # permute observations
#   indices <- sample(indices)
#   # assign a block to each observation
#   blocks <- rep(seq_len(K), length.out = n)
#   # split the permuted observations according to the block they belong to
#   folds <- split(indices, blocks)
#   names(folds) <- NULL
#   folds
# }
