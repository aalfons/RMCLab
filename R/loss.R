# ************************************
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ************************************

## pseudo-Huber loss
pseudo_huber <- function(x, delta = 1) {
  delta^2 * (sqrt(1 + (x/delta)^2) - 1)
}

## bounded absolute loss
bounded <- function(x, bound) {
  pmin(abs(x), bound)
}
