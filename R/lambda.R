# ************************************
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ************************************


#' @name lambda_grid
#'
#' @title Construct grid of values for the regularization parameter
#' 
#' @description Construct a grid of values for the regularization parameter in 
#' \code{\link{rdmc}()} or \code{\link{soft_impute}()}.
#' 
#' @details 
#' Function \code{fraction_grid()} generates a grid of values in the interval 
#' (0, 1], either on a logarithmic or linear scale, which \code{\link{rdmc}()} 
#' and \code{\link{soft_impute}()} can relate to a certain reference value 
#' computed from the data at hand. 
#' 
#' Function \code{mult_grid()} generates a multiplicative grid in which the 
#' each value is obtained by multiplying the previous value with a specified 
#' factor.
#' 
#' @param min  numeric; the smallest value of the regularization parameter.  
#' For \code{fraction_grid()}, it must be in the interval (0, 1) with the 
#' default being 0.01.  For \code{mult_grid()}, it must be larger than 0 with 
#' the default being 0.05.
#' @param max  numeric; the largest value of the regularization parameter.  It 
#' must be in the interval (\code{min}, 1] with the default being 1.
#' @param nb_lambda  a positive integer giving the number of values for the 
#' regularization parameter to be generated.
#' @param log  a logical indicating whether the grid of values should be on a 
#' logarithmic scale (defaults to \code{TRUE}).
#' @param reverse  a logical indicating whether the grid of values should be 
#' in ascending order (\code{FALSE}, the default) or in descending order 
#' (\code{TRUE}).
#' @param factor  numeric; multiplication factor larger than 1 to be used to 
#' construct the values of the regularization parameter.  That is, the second 
#' value is obtained by multiplying \code{min} by \code{factor}, with this 
#' process being iterated further.
#' 
#' @return A numeric vector of values for the regularization parameter.
#' 
#' @seealso \code{\link{rdmc}()}, \code{\link{rdmc_tune}()},
#' \code{\link{soft_impute}()}, \code{\link{soft_impute_tune}()}
#' 
#' @examples
#' fraction_grid()
#' fraction_grid(log = FALSE)
#' mult_grid(factor = 2)
#' 
#' @keywords utilities

NULL


## grid of decreasing tuning parameters in interval (0, 1], either equally 
## spaced on a logarithmic or linear scale: intended for soft_impute() as 
## fractions of softImpute::lambda0()
#' @rdname lambda_grid
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


## grid of increasing tuning parameters where the previous value is multiplied 
## with a certain factor to exponentially grow the grid: intended for rdmc()
#' @rdname lambda_grid
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
