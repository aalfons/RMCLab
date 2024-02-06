# generate data
n <- 20
p <- 5
# n <- 300
# p <- 200
nb_cat <- 5
X <- matrix(sample(1:nb_cat, n * p, replace = TRUE), nrow = n, ncol = p)
X_NA <- X
X_NA[runif(n * p) < 0.6] <- NA_real_

# other arguments
# is_NA <- is.na(X_NA)
# storage.mode(is_NA) <- "integer"
type <- "svd"
lambda <- 0.2

# loss function
pseudo_huber <- function(x, delta = 1) {
  delta^2 * (sqrt(1 + (x/delta)^2) - 1)
}
loss <- "pseudo_huber"

# R implementation
library("recommenderlab")
source("../rmc-ratings/R/rmc_ratings_funs.R")
foo <- rdmc_R(X_NA, nb_cat = nb_cat, type = type, lambda = lambda, 
              loss = "pseudo_huber")

# compare with C++ implementation from this package
bar <- rdmc(X_NA, nb_cat = nb_cat, lambda = lambda, type = type, 
            loss = "pseudo_huber")

# check if results are the same
identical(foo$L, bar$L)

# check computation time
library("microbenchmark")
microbenchmark(
  R = rdmc_R(X_NA, nb_cat = nb_cat, type = type, lambda = lambda, 
             loss = "pseudo_huber"),
  Cpp = rdmc(X_NA, nb_cat = nb_cat, lambda = lambda, type = type, 
             loss = "pseudo_huber"),
  times = 10
)
