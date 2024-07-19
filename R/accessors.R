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
get_X.rdmc_autotuned <- function(object, ...) {
  object$fit$X[[object$which_opt]]
}

#' @export
get_X.soft_impute_tuned <- function(object, discretized = FALSE, ...) {
  if (isTRUE(discretized)) object$fit$X_discretized
  else object$fit$X
}


## extract number of iterations

#' @export
get_nb_iter <- function(object, ...) UseMethod("get_nb_iter")

#' @export
get_nb_iter.rdmc_tuned <- function(object, ...) {
  object$fit$nb_iter
}

#' @export
get_nb_iter.rdmc_autotuned <- function(object, ...) {
  object$fit$nb_iter[[object$which_opt]]
}
