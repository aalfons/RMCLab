# ************************************
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ************************************


## control object for selection of tuning parameter as a fraction of lambda0()
#' @export

fraction_control <- function(start = 1, end = 0.01, nb_lambda = 10L) {
  if (start <= 0 && start > 1) {
    stop("starting fraction of tuning parameter must be in interval (0, 1]")
  }
  if (end <= 0 && end >= start) {
    stop("ending fraction of tuning parameter must be in interval (0, start)")
  }
  if (nb_lambda <= 0) {
    stop("number of tuning parameter values must be greater than 0")
  }
  out <- list(start = start, end = end, nb_lambda = nb_lambda)
  class(out) <- "fraction_control"
  out
}


# control object for selection
rel_lambda_control <- function(smallest = 0.05, factor = 1.5, nb_lambda = 10L) {
  if (smallest <= 0) {
    stop("smallest value of tuning parameter must be greater than 0")
  }
  if (factor <= 1) {
    stop("multiplication factor must be greater than 1")
  }
  if (nb_lambda <= 0) {
    stop("number of tuning parameter values must be greater than 0")
  }
  out <- list(smallest = smallest, factor = factor, nb_lambda = nb_lambda)
  class(out) <- "rel_lambda_control"
  out
}

## obtain grid of tuning parameter values

#' @noRd
get_grid <- function(object, ...) UseMethod("get_grid")

#' @noRd
get_grid.fraction_control <- function(object, ...) {
  seq_log <- seq(from = log(object$start), to = log(object$end),
                 length.out = object$nb_lambda)
  exp(seq_log)
}

#' @noRd
get_grid.rel_lambda_control <- function(object, ...) {
  seq_exponent <- seq(from = 0L, to = object$nb_lambda - 1L)
  object$smallest * object$factor^seq_exponent
}
