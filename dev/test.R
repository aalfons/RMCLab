library("rdmc")

# set seed of random number generator for reproducibility
set.seed(1234)

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
lambda <- c(0.01, 0.2, 0.4, 0.6, 0.8)

# single lambda
fit <- rdmc(X_NA, values = 1:nb_cat, lambda = lambda,
            type = type, loss = "pseudo_huber")
fit$lambda
fit$d_max
fit$nb_iter
fit$L
fit$X

foo <- rdmc(X_NA, values = 1:nb_cat, lambda = lambda * fit$d_max, 
            relative = FALSE, type = type, loss = "pseudo_huber")
foo$lambda
foo$d_max
identical(fit$X, foo$X)

# pseudo-Huber loss
fit_pseudo_huber <- rdmc_tune(X_NA, values = 1:nb_cat, lambda = lambda, 
                              type = type, loss = "pseudo_huber")
fit_pseudo_huber$lambda_opt
fit_pseudo_huber$fit$X

# absolute loss
fit_absolute <- rdmc_tune(X_NA, values = 1:nb_cat, lambda = lambda, type = type, 
                          loss = "absolute")
fit_absolute$lambda_opt
fit_absolute$fit$X

# truncated absolute loss
fit_truncated <- rdmc_tune(X_NA, values = 1:nb_cat, lambda = lambda, type = type, 
                           loss = "truncated")
fit_truncated$lambda_opt
fit_truncated$fit$X
