# ************************************
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ************************************


## extract imputed data matrix

#' @export
get_X <- function(object, ...) UseMethod("get_X")

#' @export
get_X.rdmc_tuned <- function(object, ...) {
  object$fit$X
}

#' @export
get_X.soft_impute_tuned <- function(object, discretized = FALSE, ...) {
  if (isTRUE(discretized)) object$fit$X_discretized
  else object$fit$X
}

#' @export
get_X.median_impute <- function(object, discretized = FALSE, ...) {
  if (isTRUE(discretized)) object$X_discretized
  else object$X
}

#' @export
get_X.mode_impute <- function(object, ...) {
  object$X
}


## extract value of the tuning parameter

#' @export
get_lambda <- function(object, ...) UseMethod("get_lambda")

#' @export
get_lambda.rdmc <- function(object, relative = TRUE, ...) {
  lambda <- object$lambda
  if (isTRUE(relative)) lambda 
  else lambda * object$d_max
}

#' @export
get_lambda.rdmc_tuned <- function(object, ...) {
  get_lambda(object$fit, ...)
}

#' @export
get_lambda.soft_impute <- function(object, relative = TRUE, ...) {
  lambda <- object$lambda
  if (isTRUE(relative)) lambda 
  else lambda * object$lambda0
}

#' @export
get_lambda.soft_impute_tuned <- function(object, ...) {
  get_lambda(object$fit, ...)
}


## extract number of iterations

#' @export
get_nb_iter <- function(object, ...) UseMethod("get_nb_iter")

#' @export
get_nb_iter.rdmc_tuned <- function(object, ...) {
  object$fit$nb_iter
}
