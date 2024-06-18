library("rdmc")

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
lambda <- c(0.2, 0.4, 0.6, 0.8, 1)

# pseudo-Huber loss
fit_pseudo_huber <- rdmc_tune(X_NA, values = 1:nb_cat, lambda = lambda, 
                              type = type, loss = "pseudo_huber")
fit_pseudo_huber$lambda_opt
fit_pseudo_huber$final$X

# absolute loss
fit_absolute <- rdmc_tune(X_NA, values = 1:nb_cat, lambda = lambda, type = type, 
                          loss = "absolute")
fit_absolute$lambda_opt
fit_absolute$final$X

# bounded absolute loss
fit_bounded <- rdmc_tune(X_NA, values = 1:nb_cat, lambda = lambda, type = type, 
                         loss = "bounded")
fit_bounded$lambda_opt
fit_bounded$final$X
