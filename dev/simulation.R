# ***************************************
# Authors: Andreas Alfons
#          Erasmus Universiteit Rotterdam
#
#          Aurore Archimbaud
#          Toulouse Business School
# ***************************************


# control paramaters for data generation
n <- 300      # number of observations
p <- 200      # number of variables 
J <- 50       # rank
sigma <- 0.2  # standard deviation of noise

# control parameters for discretization
nb_cat <- 5                                             # number of categories
values <- 1:nb_cat                                      # rating-scale values
breaks <- c(-Inf, (values[-1]+values[-nb_cat])/2, Inf)  # cut points

# control parameters for missing values
NA_fractions <- c(0.6, 0.9)

# other control parameters
R <- 50           # number of simulation runs
seed <- 20240207  # seed of the random number generator

# control parameters for methods
svd_tol <- 1e-4
delta <- 1.05 
mu <- 0.1
conv_tol <- 1e-4
max_iter <- 100


# evaluation measures
MAE <- function(true, predicted) mean(abs(true - predicted))
RMSE <- function(true, predicted) sqrt(mean((true - predicted)^2))


# for parallel computing
n_cores <- 1
# if (.Platform$OS.type == "windows") {
#   n_cores <- 1              # use only one CPU core on Windows
# } else {
#   n_cores <- 2              # number of CPU cores to be used
#   RNGkind("L'Ecuyer-CMRG")  # use parallel random number streams
# }


# print message that simulation is starting
cat(paste(Sys.time(), ": starting ...\n"))

# start simulation
set.seed(seed)
results_list <- parallel::mclapply(seq_len(R), function(r) {
  
  # print simulation run
  cat(paste(Sys.time(), sprintf(":   run = %d\n", r)))
  
  # generate data from Gaussian factor model with Gaussian noise
  M_continuous <- 
    matrix(rnorm(n*J), nrow = n, ncol = J) %*% 
    matrix(rnorm(J*p), nrow = J, ncol = p) +
    sigma * matrix(rnorm(n*p), nrow = n, ncol = p)
  
  # discretize the generated data
  M_shifted <- scale(M_continuous, center = TRUE, scale = TRUE) + (1+nb_cat)/2
  X <- apply(M_shifted, 2, function(x) as.numeric(cut(x, breaks = breaks)))
  
  # loop over missingness fractions
  results_NA_fractions <- lapply(NA_fractions, function(NA_frac) {
    
    # replace observations with missing values
    set_NA <- sample.int(n*p, NA_frac * n*p)
    X_NA <- X
    X_NA[set_NA] <- NA_real_
    
    # grid of values of for regularization parameter
    # (note that we soft-threshold the singluar values with lambda / mu, while 
    # softImpute soft-thresholds with lambda, so we need to correct the estimate 
    # of the smallest lambda that sets all singular values to 0)
    lambda0 <- softImpute::lambda0(X_NA) * mu
    lambda_grid <- exp(seq(from = log(1), to = log(lambda0), length.out = 10))
    
    # rdmc() without warm starts
    df_cold_starts <- tryCatch({
      out_list <- lapply(lambda_grid, function(lambda) {
        rdmc::rdmc(X_NA, nb_cat = nb_cat, lambda = lambda, 
                   svd_tol = svd_tol, delta = delta, mu = mu, 
                   conv_tol = conv_tol, max_iter = max_iter)
      })
      mae <- sapply(out_list, function(out) MAE(X[set_NA], out$L[set_NA]))
      rmse <- sapply(out_list, function(out) RMSE(X[set_NA], out$L[set_NA]))
      data.frame(Run = r, n = n, p = p, J = J, nb_cat = nb_cat, 
                 pct_missing = NA_frac, Method = "RDMC_cold_starts",
                 lambda = seq_along(lambda_grid), MAE = mae / (nb_cat-1), 
                 RMSE = rmse / (nb_cat-1))
    }, error = function(e) NULL)
    
    # rdmc() with warm starts
    df_warm_starts <- tryCatch({
      out_list <- rdmc::rdmc(X_NA, nb_cat = nb_cat, lambda = lambda_grid, 
                             svd_tol = svd_tol, delta = delta, mu = mu, 
                             conv_tol = conv_tol, max_iter = max_iter)
      mae <- sapply(out_list, function(out) MAE(X[set_NA], out$L[set_NA]))
      rmse <- sapply(out_list, function(out) RMSE(X[set_NA], out$L[set_NA]))
      data.frame(Run = r, n = n, p = p, J = J, nb_cat = nb_cat, 
                 pct_missing = NA_frac, Method = "RDMC_warm_starts",
                 lambda = seq_along(lambda_grid), MAE = mae / (nb_cat-1), 
                 RMSE = rmse / (nb_cat-1))
    }, error = function(e) NULL)
    
    # combine results from current simulation run into data frame
    rbind(df_cold_starts, df_warm_starts)
    
  })
  
  # combine results for different missingness fractions into data frame
  do.call(rbind, results_NA_fractions)
  
}, mc.cores = n_cores)


# combine results into data frame
results <- do.call(rbind, results_list)

# save results to file
file_results <- "dev/results_n=%d_p=%d_J=%d_R=%d.RData"
save(results, seed, file = sprintf(file_results, n, p, J, R))

# print message that simulation is done
cat(paste(Sys.time(), ": finished.\n"))


# load packages for plotting results
library("dplyr")
library("ggplot2")

# prepare data frame
df_plot <- results %>%
  mutate(lambda = as.factor(lambda))

# generate plot of MAE
ggplot() +
  geom_boxplot(aes(x = lambda, y = MAE, fill = Method),
               data = df_plot) +
  facet_grid(. ~ pct_missing) +
  theme(legend.position = "top")

# generate plot of RMSE
ggplot() +
  geom_boxplot(aes(x = lambda, y = RMSE, fill = Method),
               data = df_plot) +
  facet_grid(. ~ pct_missing) +
  theme(legend.position = "top")
