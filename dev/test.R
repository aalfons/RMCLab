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
lambda <- 1

# loss function
pseudo_huber <- function(x, delta = 1) {
  delta^2 * (sqrt(1 + (x/delta)^2) - 1)
}

# pseudo-Huber loss
fit_pseudo_huber <- rdmc(X_NA, values = 1:nb_cat, lambda = lambda, type = type, 
                         loss = "pseudo_huber")
fit_pseudo_huber$X

# absolute loss
fit_absolute <- rdmc(X_NA, values = 1:nb_cat, lambda = lambda, type = type, 
                     loss = "absolute")
fit_absolute$X

# bounded absolute loss
fit_bounded <- rdmc(X_NA, values = 1:nb_cat, lambda = lambda, type = type, 
                     loss = "bounded")
fit_bounded$X
