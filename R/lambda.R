# ************************************
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ************************************


## grid of increasing tuning parameters where the previous value is multiplied 
## with a certain factor to exponentially grow the grid: intended for rdmc()
#' @export
mult_grid <- function(smallest = 0.05, factor = 1.5, nb_lambda = 10L) {
  if (smallest <= 0) {
    stop("smallest value of tuning parameter must be greater than 0")
  }
  if (factor <= 1) {
    stop("multiplication factor must be greater than 1")
  }
  if (nb_lambda <= 0) {
    stop("number of tuning parameter values must be greater than 0")
  }
  seq_exponent <- seq(from = 0L, to = nb_lambda - 1L)
  smallest * factor^seq_exponent
}


## grid of decreasing tuning parameters in interval (0, 1], either equally 
## spaced on a logarithmic or linear scale: intended for soft_impute() as 
## fractions of softImpute::lambda0()
#' @export
fraction_grid <- function(start = 1, end = 0.01, nb_lambda = 10L, log = TRUE) {
  if (start <= 0 && start > 1) {
    stop("starting fraction of tuning parameter must be in interval (0, 1]")
  }
  if (end <= 0 && end >= start) {
    stop("ending fraction of tuning parameter must be in interval (0, start)")
  }
  if (nb_lambda <= 0) {
    stop("number of tuning parameter values must be greater than 0")
  }
  if (isTRUE(log)) {
    seq_log <- seq(from = log(start), to = log(end), length.out = nb_lambda)
    exp(seq_log)
  } else {
    seq(from = start, to = end, length.out = nb_lambda)
  }
}
