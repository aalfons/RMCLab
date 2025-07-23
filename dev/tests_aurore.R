# Installation ----
#devtools::install_github("aalfons/rdmc", ref = "dev")
#library(rdmc)


# Title: should we say imputation of rating scale data either in the title or the description?

# Import data  ----
load("data/MovieLens_subset.Rdata")
# add data? and description? name of data? X_NA?
MovieLensToy

# rda or Rdata?
#usethis::use_data("data/HTP2.rda")

# Imputation -----
## median -----

fit_median <- median_impute(X_NA, values = 1:5)

# I added a message : "You requested the discretized imputed matrix, but this step must be explicitly enabled in the median_impute() function. 
# Please set the parameter discretize = TRUE to perform this step."
# By default it is true
X_median <- get_completed(fit_median)
X_median <- get_imputed(fit_median)

fit_median <- median_impute(X_NA, discretize = FALSE, values = 1:5)
fit_median_disc <- median_impute(X_NA, discretize = TRUE, values = 1:5)
X_median <- get_completed(fit_median)


## mode -----
fit_mode <- mode_impute(X_NA, values = 1:5)
X_mode <- get_completed(fit_mode)


## rdmc -----
fit_rdmc <- rdmc(X_NA)
fit_rdmc <- rdmc(X_NA, loss = "pseudo_huber")
X_rdmc <- get_completed(fit_rdmc)
X_rmdc  <- get_completed(fit_rdmc, which = 1)
get_nb_iter(fit_rdmc_tune)
get_lambda(fit_rdmc_tune)
## rdmc_tune -----
# Automatically choose the best regularization parameter

# By default: selection via repeated holdout validation.
fit_rdmc_tuned <- rdmc_tune(X_NA, loss = "truncated")
X_rdmc_tuned <- get_completed(fit_rdmc_tuned)
get_nb_iter(fit_rdmc_tuned)
get_lambda(fit_rdmc_tuned)

# Warning: nothing by default. I think it is because you added the previous bounded and truncated one. 

# Change lambda grid and splits
lambda <- fraction_grid(min = 0.01, max = 1, nb_lambda = 7, log = TRUE)
splits <- holdout_control(pct = 0.1, R = 5)
fit_rdmc_tuned <- rdmc_tune(X_NA, loss = "absolute", values =  1:5, 
                           lambda = lambda, splits = splits,
                           conv_tol =  1e-4, max_iter = 10)
X_rdmc_tuned <- get_completed(fit_rdmc_tuned)
get_nb_iter(fit_rdmc_tuned)
get_lambda(fit_rdmc_tuned)



## soft_impute -----
fit_soft_impute <- soft_impute(X_NA)
X_soft_impute <- get_completed(fit_soft_impute, which = 1)

# by default no discretization


## soft_impute_tune -----
fit_soft_impute_tuned <- soft_impute_tune(X_NA, discretize = TRUE)

X_soft_impute_tuned <- get_completed(fit_soft_impute_tuned, discretized = FALSE)
X_soft_impute_tuned_discretized <- get_completed(fit_soft_impute_tuned, 
                                                discretized = TRUE)
get_lambda(fit_soft_impute_tuned)



lambda <- fraction_grid(min = 0.01, max = 1, nb_lambda = 7, log = TRUE)
splits <- holdout_control(pct = 0.1, R = 5)
fit_soft_impute_tune <- soft_impute_tune(X_NA,
                                         lambda = lambda, splits = splits,
                                         thresh = 1e-04,
                                         maxit = 10)
X_soft_impute_tuned <- get_completed(fit_soft_impute_tuned, discretized = FALSE)
X_soft_impute_tuned <- get_completed(fit_soft_impute_tuned, discretized = TRUE)
get_lambda(fit_soft_impute_tuned)
# by default discretization

# Create splits of observed data cells for hyperparameter tuning
create_splits(indices, control)
holdout(indices, pct = 0.1, R = 10L)
cv_folds(indices, K = 5L)

#Extract the optimal value of the regularization parameter
get_lambda 
# not clear

#Construct grid of values for the regularization parameter
lambda_grid 
fraction_grid()
fraction_grid(log = FALSE)
mult_grid(factor = 2)
# maybe add a sentence for the description

# New repo with new figures as well.
