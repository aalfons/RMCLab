# ************************************
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ************************************

## pseudo-Huber loss
pseudo_huber <- function(x, const = 1) {
  const^2 * (sqrt(1 + (x/const)^2) - 1)
}

## truncated absolute loss
truncated <- function(x, const) {
  pmin(abs(x), const)
}
