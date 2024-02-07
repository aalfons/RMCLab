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
type <- "svd"
lambda <- 0.2
max_iter = 100

# loss function
pseudo_huber <- function(x, delta = 1) {
  delta^2 * (sqrt(1 + (x/delta)^2) - 1)
}
loss <- "pseudo_huber"

# R implementation
library("recommenderlab")
source("../rmc-ratings/R/rmc_ratings_funs.R")
foo <- rdmc_R(X_NA, nb_cat = nb_cat, type = type, lambda = lambda, 
              loss = loss, max_iter = max_iter)

# compare with C++ implementation from this package
library("rdmc")
bar <- rdmc(X_NA, nb_cat = nb_cat, lambda = lambda, type = type, 
            loss = loss, max_iter = max_iter)

# check if results are the same
identical(foo$L, bar$L)

# check number of iterations
bar$nb_iter

# check computation time
library("microbenchmark")
microbenchmark(
  R = rdmc_R(X_NA, nb_cat = nb_cat, type = type, lambda = lambda, 
             loss = loss, max_iter = max_iter),
  Cpp = rdmc(X_NA, nb_cat = nb_cat, lambda = lambda, type = type, 
             loss = loss, max_iter = max_iter),
  times = 10
)


# grid of values for the regularization parameter
lambdas <- seq(from = 0, to = 1, by = 0.1)

# R implementation without warm starts
res_R <- lapply(lambdas, function(lambda) {
  rdmc_R(X_NA, nb_cat = nb_cat, type = type, lambda = lambda, 
         loss = loss,  max_iter = max_iter)
})

# C++ implementation without warm starts
res_cold_starts <- lapply(lambdas, function(lambda) {
  rdmc(X_NA, nb_cat = nb_cat, lambda = lambda, type = type, 
       loss = loss, max_iter = max_iter)
})

# C++ implementation with warm starts
res_warm_starts <- rdmc(X_NA, nb_cat = nb_cat, lambda = lambdas, type = type, 
                        loss = loss, max_iter = max_iter)

# check computation time
library("microbenchmark")
microbenchmark(
  R_full = lapply(lambdas, function(lambda) {
    rdmc_R(X_NA, nb_cat = nb_cat, type = type, lambda = lambda, 
           loss = loss, max_iter = max_iter)
  }),
  R_10_iter = lapply(lambdas, function(lambda) {
    rdmc_R(X_NA, nb_cat = nb_cat, type = type, lambda = lambda, 
           loss = loss, max_iter = 10)
  }),
  Cpp_cold_starts = lapply(lambdas, function(lambda) {
    rdmc(X_NA, nb_cat = nb_cat, lambda = lambda, type = type, 
         loss = loss, max_iter = max_iter)
  }),
  Cpp_warm_starts = rdmc(X_NA, nb_cat = nb_cat, lambda = lambdas, type = type, 
                         loss = loss, max_iter = max_iter),
  times = 10
)
