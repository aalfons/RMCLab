# ************************************
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ************************************


## grid of increasing tuning parameters where the previous value is multiplied 
## with a certain factor to exponentially grow the grid: intended for rdmc()
#' @export
mult_grid <- function(min = 0.05, factor = 1.5, nb_lambda = 10L) {
  if (min <= 0) {
    stop("smallest value of tuning parameter must be greater than 0")
  }
  if (factor <= 1) {
    stop("multiplication factor must be greater than 1")
  }
  if (nb_lambda <= 0) {
    stop("number of tuning parameter values must be greater than 0")
  }
  seq_exponent <- seq(from = 0L, to = nb_lambda - 1L)
  min * factor^seq_exponent
}


## grid of decreasing tuning parameters in interval (0, 1], either equally 
## spaced on a logarithmic or linear scale: intended for soft_impute() as 
## fractions of softImpute::lambda0()
#' @export
fraction_grid <- function(min = 0.01, max = 1, nb_lambda = 10L, 
                          log = TRUE, reverse = FALSE) {
  if (min <= 0 || min >= 1) {
    stop("smallest fraction of tuning parameter must be in interval (0, 1)")
  }
  if (max <= min || max > 1) {
    stop("largest fraction of tuning parameter must be in interval (min, 1]")
  }
  if (nb_lambda <= 0) {
    stop("number of tuning parameter values must be greater than 0")
  }
  if (isTRUE(log)) {
    seq_log <- seq(from = log(min), to = log(max), length.out = nb_lambda)
    grid <- exp(seq_log)
  } else {
    grid <- seq(from = min, to = max, length.out = nb_lambda)
  }
  if (isTRUE(reverse)) rev(grid)
  else grid
}
