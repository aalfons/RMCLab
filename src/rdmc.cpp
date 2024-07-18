/*
 * Author: Andreas Alfons
 *         Erasmus Universiteit Rotterdam
 * 
 * based on R code by Aurore Archimbaud
 */

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>


// For now, there are multiple copies of the function where the only difference
// is how the loss function is computed in step 2 (specifically, update L for 
// cells with observed values in X). This is because the loss function is 
// computed for every observed cell in every iteration, so we should avoid if 
// statements inside those loops. We should look into object-oriented 
// implementations using class templates, which could be an efficient way of 
// handling the loss function inside the loops. That is, we should have 
// different classes for the different loss function, with a method that 
// computes the loss for a given cell.


// -----------------
// pseudo-Huber loss
// -----------------

// loss function
double pseudo_huber(const double& x, const double& delta) {
  return std::pow(delta, 2.0) * (sqrt(1 + std::pow(x/delta, 2.0)) - 1);
}

// workhorse function for a single value of the regularization parameter lambda
void rdmc_pseudo_huber(const arma::mat& X, const arma::uword& n, 
                       const arma::uword& p, const arma::uvec& ind_NA, 
                       const arma::uvec& ind_not_NA, const arma::mat& values, 
                       const double& lambda, const std::string& type, 
                       const double& svd_tol, const double& loss_tuning, 
                       const double& delta, double mu, const double& conv_tol, 
                       const arma::uword& max_iter,
                       // output to be returned through arguments
                       arma::mat& L, arma::mat& Z, arma::mat& Theta,
                       double& objective, bool& converged, 
                       arma::uword& nb_iter) {

  // initializations
  objective = R_PosInf;
  converged = false;
  nb_iter = 0;
  
  // iterate update steps of Z, L, and Theta
  arma::uword rank, i, j, k, replace_i, nb_values = values.n_rows, which_min;
  arma::mat U, V, L_minus_Z;
  arma::vec d;
  double nuclear_norm, tmp, objective_step2, objective_step2_min,
  loss_norm, loss, loss_min, previous_objective, change;
  while (!converged && nb_iter < max_iter) {
    
    // step 1: update Z keeping L fixed
    // soft-thresholding of L + 1/mu * Theta
    
    // compute SVD
    // TODO: add option to use SVD-ALS or softImpute-ALS instead of the SVD
    arma::svd(U, d, V, L + Theta/mu);
    // adjust singular values for our parametrization
    d -= lambda / mu;
    // compute rank of soft-thresholded singular values
    nuclear_norm = 0;
    rank = 0;
    for (j = 0; j < d.n_elem; j++) {
      if (d(j) > svd_tol) {
        nuclear_norm += d(j);
        rank++;
      }
    }
    // initialize Z with zeros
    Z.zeros();
    // if we have a nonzero soft-thresholded singular value, update Z from
    // the soft-thresholded SVD
    if (rank > 0) {
      // efficient implementation of matrix multiplications without copying parts
      for (i = 0; i < n; i++) {
        for (j = 0; j < p; j++) {
          // add contributions from the different rows of U and columns of V'
          for (k = 0; k < rank; k++) Z(i, j) += U(i, k) * d(k) * V(j, k);
        }
      }
    }
    
    // step 2: update L keeping Z fixed
    // separable problem for missing and observed values in X
    
    // update L for cells with missing values in X
    for (i = 0; i < ind_NA.n_elem; i++) {
      // index of current cell to be updated
      replace_i = ind_NA(i);
      // column of current cell to be updated
      // TODO: this should be avoided by passing down indices (i, j) 
      // for which cells are missing
      j = replace_i / n;
      // save some computation time in loop
      tmp = -Z(replace_i) + Theta(replace_i)/mu;
      // initialize the minimum
      which_min = 0;
      objective_step2_min = R_PosInf;
      // loop over the different values and choose the one that minimizes the
      // objective function
      for (k = 0; k < nb_values; k++) {
        // The paper says to take the argmin of the squared expression. But
        // this is equivalent to taking the argmin of the absolute value,
        // which is faster to compute.
        objective_step2 = std::abs(values(k, j) + tmp);
        if (objective_step2 < objective_step2_min) {
          which_min = k;
          objective_step2_min = objective_step2;
        }
      }
      // update the element of L with the argmin of the objective function
      L(replace_i) = values(which_min, j);
    }
    
    // update L for cells with observed values in X
    loss_norm = 0;
    for (i = 0; i < ind_not_NA.n_elem; i++) {
      // index of current cell to be updated
      replace_i = ind_not_NA(i);
      // column of current cell to be updated
      // TODO: this should be avoided by passing down indices (i, j) 
      // for which cells are missing
      j = replace_i / n;
      // save some computation time in loop
      tmp = -Z(replace_i) + Theta(replace_i)/mu;
      // initialize the minimum
      which_min = 0;
      objective_step2_min = R_PosInf;
      // loop over the different values and choose the one that minimizes the
      // objective function
      for (k = 0; k < nb_values; k++) {
        loss = pseudo_huber(values(k, j) - X(replace_i), loss_tuning);
        objective_step2 = loss + mu * std::pow(values(k, j) + tmp, 2.0)/2.0;
        if (objective_step2 < objective_step2_min) {
          which_min = k;
          loss_min = loss;
          objective_step2_min = objective_step2;
        }
      }
      // update the element of L with the argmin of the objective function
      L(replace_i) = values(which_min, j);
      // update the norm given by loss function
      loss_norm += loss_min;
    }
    
    // step 3: update Lagrange multiplier Theta and parameter mu
    L_minus_Z = L - Z;
    Theta = Theta + mu * L_minus_Z;
    mu = delta * mu;
    
    // update iteration counter
    nb_iter++;
    // update objective function for convergence criterion
    previous_objective = objective;
    objective = loss_norm + lambda * nuclear_norm + arma::dot(Theta, L_minus_Z) +
      mu * arma::norm(L_minus_Z, "fro") / 2.0;
    // compute relative change and check convergence
    if (nb_iter > 1) {
      // we can't compute relative change in the first iteration since the
      // objective function is initialized with infinity
      change = std::abs((objective - previous_objective) / previous_objective);
      converged = change < conv_tol;
    }
    
  }
  
}


// workhorse function for a single value of the regularization parameter lambda
void rdmc_pseudo_huber_rank_max(const arma::mat& X, const arma::uword& n, 
                                const arma::uword& p, const arma::uvec& ind_NA, 
                                const arma::uvec& ind_not_NA, const arma::mat& values, 
                                const double& lambda, const arma::uword& rank_max, 
                                const std::string& type, const double& svd_tol, 
                                const double& loss_tuning, const double& delta, 
                                double mu, const double& conv_tol, 
                                const arma::uword& max_iter,
                                // output to be returned through arguments
                                arma::mat& L, arma::mat& Z, arma::mat& Theta,
                                double& objective, bool& converged, 
                                arma::uword& nb_iter) {
  
  // initializations
  objective = R_PosInf;
  converged = false;
  nb_iter = 0;
  
  // iterate update steps of Z, L, and Theta
  arma::uword rank, i, j, k, replace_i, nb_values = values.n_rows, which_min;
  arma::mat U, V, L_minus_Z;
  arma::vec d;
  double nuclear_norm, tmp, objective_step2, objective_step2_min,
  loss_norm, loss, loss_min, previous_objective, change;
  while (!converged && nb_iter < max_iter) {
    
    // step 1: update Z keeping L fixed
    // soft-thresholding of L + 1/mu * Theta
    
    // compute SVD
    // TODO: add option to use SVD-ALS or softImpute-ALS instead of the SVD
    arma::svd(U, d, V, L + Theta/mu);
    // adjust singular values for our parametrization
    d -= lambda / mu;
    // compute rank of soft-thresholded singular values
    nuclear_norm = 0;
    rank = 0;
    // for (j = 0; j < d.n_elem; j++) {
    //   if (d(j) > svd_tol) {
    //     nuclear_norm += d(j);
    //     rank++;
    //   }
    // }
    j = 0;
    while ((j < d.n_elem) & (rank < rank_max)) {
      if (d(j) > svd_tol) {
        nuclear_norm += d(j);
        rank++;
      }
      j++;
    }
    // initialize Z with zeros
    Z.zeros();
    // if we have a nonzero soft-thresholded singular value, update Z from
    // the soft-thresholded SVD
    if (rank > 0) {
      // efficient implementation of matrix multiplications without copying parts
      for (i = 0; i < n; i++) {
        for (j = 0; j < p; j++) {
          // add contributions from the different rows of U and columns of V'
          for (k = 0; k < rank; k++) Z(i, j) += U(i, k) * d(k) * V(j, k);
        }
      }
    }
    
    // step 2: update L keeping Z fixed
    // separable problem for missing and observed values in X
    
    // update L for cells with missing values in X
    for (i = 0; i < ind_NA.n_elem; i++) {
      // index of current cell to be updated
      replace_i = ind_NA(i);
      // column of current cell to be updated
      // TODO: this should be avoided by passing down indices (i, j) 
      // for which cells are missing
      j = replace_i / n;
      // save some computation time in loop
      tmp = -Z(replace_i) + Theta(replace_i)/mu;
      // initialize the minimum
      which_min = 0;
      objective_step2_min = R_PosInf;
      // loop over the different values and choose the one that minimizes the
      // objective function
      for (k = 0; k < nb_values; k++) {
        // The paper says to take the argmin of the squared expression. But
        // this is equivalent to taking the argmin of the absolute value,
        // which is faster to compute.
        objective_step2 = std::abs(values(k, j) + tmp);
        if (objective_step2 < objective_step2_min) {
          which_min = k;
          objective_step2_min = objective_step2;
        }
      }
      // update the element of L with the argmin of the objective function
      L(replace_i) = values(which_min, j);
    }
    
    // update L for cells with observed values in X
    loss_norm = 0;
    for (i = 0; i < ind_not_NA.n_elem; i++) {
      // index of current cell to be updated
      replace_i = ind_not_NA(i);
      // column of current cell to be updated
      // TODO: this should be avoided by passing down indices (i, j) 
      // for which cells are missing
      j = replace_i / n;
      // save some computation time in loop
      tmp = -Z(replace_i) + Theta(replace_i)/mu;
      // initialize the minimum
      which_min = 0;
      objective_step2_min = R_PosInf;
      // loop over the different values and choose the one that minimizes the
      // objective function
      for (k = 0; k < nb_values; k++) {
        loss = pseudo_huber(values(k, j) - X(replace_i), loss_tuning);
        objective_step2 = loss + mu * std::pow(values(k, j) + tmp, 2.0)/2.0;
        if (objective_step2 < objective_step2_min) {
          which_min = k;
          loss_min = loss;
          objective_step2_min = objective_step2;
        }
      }
      // update the element of L with the argmin of the objective function
      L(replace_i) = values(which_min, j);
      // update the norm given by loss function
      loss_norm += loss_min;
    }
    
    // step 3: update Lagrange multiplier Theta and parameter mu
    L_minus_Z = L - Z;
    Theta = Theta + mu * L_minus_Z;
    mu = delta * mu;
    
    // update iteration counter
    nb_iter++;
    // update objective function for convergence criterion
    previous_objective = objective;
    objective = loss_norm + lambda * nuclear_norm + arma::dot(Theta, L_minus_Z) +
      mu * arma::norm(L_minus_Z, "fro") / 2.0;
    // compute relative change and check convergence
    if (nb_iter > 1) {
      // we can't compute relative change in the first iteration since the
      // objective function is initialized with infinity
      change = std::abs((objective - previous_objective) / previous_objective);
      converged = change < conv_tol;
    }
    
  }
  
}


// -------------
// absolute loss
// -------------

// workhorse function for a single value of the regularization parameter lambda
void rdmc_absolute(const arma::mat& X, const arma::uword& n,
                  const arma::uword& p, const arma::uvec& ind_NA,
                  const arma::uvec& ind_not_NA, const arma::mat& values,
                  const double& lambda, const std::string& type,
                  const double& svd_tol, const double& delta, double mu,
                  const double& conv_tol, const arma::uword& max_iter,
                  // output to be returned through arguments
                   arma::mat& L, arma::mat& Z, arma::mat& Theta,
                   double& objective, bool& converged,
                   arma::uword& nb_iter) {

  // initializations
  objective = R_PosInf;
  converged = false;
  nb_iter = 0;

  // iterate update steps of Z, L, and Theta
  arma::uword rank, i, j, k, replace_i, nb_values = values.n_rows, which_min;
  arma::mat U, V, L_minus_Z;
  arma::vec d;
  double nuclear_norm, tmp, objective_step2, objective_step2_min,
  loss_norm, loss, loss_min, previous_objective, change;
  while (!converged && nb_iter < max_iter) {

    // step 1: update Z keeping L fixed
    // soft-thresholding of L + 1/mu * Theta

    // compute SVD
    // TODO: add option to use SVD-ALS or softImpute-ALS instead of the SVD
    arma::svd(U, d, V, L + Theta/mu);
    // adjust singular values for our parametrization
    d -= lambda / mu;
    // compute rank of soft-thresholded singular values
    nuclear_norm = 0;
    rank = 0;
    for (j = 0; j < d.n_elem; j++) {
      if (d(j) > svd_tol) {
        nuclear_norm += d(j);
        rank++;
      }
    }
    // initialize Z with zeros
    Z.zeros();
    // if we have a nonzero soft-thresholded singular value, update Z from
    // the soft-thresholded SVD
    if (rank > 0) {
      // efficient implementation of matrix multiplications without copying parts
      for (i = 0; i < n; i++) {
        for (j = 0; j < p; j++) {
          // add contributions from the different rows of U and columns of V'
          for (k = 0; k < rank; k++) Z(i, j) += U(i, k) * d(k) * V(j, k);
        }
      }
    }

    // step 2: update L keeping Z fixed
    // separable problem for missing and observed values in X

    // update L for cells with missing values in X
    for (i = 0; i < ind_NA.n_elem; i++) {
      // index of current cell to be updated
      replace_i = ind_NA(i);
      // column of current cell to be updated
      // TODO: this should be avoided by passing down indices (i, j) 
      // for which cells are missing
      j = replace_i / n;
      // save some computation time in loop
      tmp = -Z(replace_i) + Theta(replace_i)/mu;
      // initialize the minimum
      which_min = 0;
      objective_step2_min = R_PosInf;
      // loop over the different values and choose the one that minimizes the
      // objective function
      for (k = 0; k < nb_values; k++) {
        // The paper says to take the argmin of the squared expression. But
        // this is equivalent to taking the argmin of the absolute value,
        // which is faster to compute.
        objective_step2 = std::abs(values(k, j) + tmp);
        if (objective_step2 < objective_step2_min) {
          which_min = k;
          objective_step2_min = objective_step2;
        }
      }
      // update the element of L with the argmin of the objective function
      L(replace_i) = values(which_min, j);
    }

    // update L for cells with observed values in X
    loss_norm = 0;
    for (i = 0; i < ind_not_NA.n_elem; i++) {
      // index of current cell to be updated
      replace_i = ind_not_NA(i);
      // column of current cell to be updated
      // TODO: this should be avoided by passing down indices (i, j) 
      // for which cells are missing
      j = replace_i / n;
      // save some computation time in loop
      tmp = -Z(replace_i) + Theta(replace_i)/mu;
      // initialize the minimum
      which_min = 0;
      objective_step2_min = R_PosInf;
      // loop over the different values and choose the one that minimizes the
      // objective function
      for (k = 0; k < nb_values; k++) {
        loss = std::abs(values(k, j) - X(replace_i));
        objective_step2 = loss + mu * std::pow(values(k, j) + tmp, 2.0)/2.0;
        if (objective_step2 < objective_step2_min) {
          which_min = k;
          loss_min = loss;
          objective_step2_min = objective_step2;
        }
      }
      // update the element of L with the argmin of the objective function
      L(replace_i) = values(which_min, j);
      // update the norm given by loss function
      loss_norm += loss_min;
    }

    // step 3: update Lagrange multiplier Theta and parameter mu
    L_minus_Z = L - Z;
    Theta = Theta + mu * L_minus_Z;
    mu = delta * mu;

    // update iteration counter
    nb_iter++;
    // update objective function for convergence criterion
    previous_objective = objective;
    objective = loss_norm + lambda * nuclear_norm + arma::dot(Theta, L_minus_Z) +
      mu * arma::norm(L_minus_Z, "fro") / 2.0;
    // compute relative change and check convergence
    if (nb_iter > 1) {
      // we can't compute relative change in the first iteration since the
      // objective function is initialized with infinity
      change = std::abs((objective - previous_objective) / previous_objective);
      converged = change < conv_tol;
    }

  }

}


// ---------------------
// bounded absolute loss
// ---------------------

// loss function
double bounded(const double& x, const double& bound) {
  return std::min(std::abs(x), bound);
}

// workhorse function for a single value of the regularization parameter lambda
void rdmc_bounded(const arma::mat& X, const arma::uword& n,
                  const arma::uword& p, const arma::uvec& ind_NA,
                  const arma::uvec& ind_not_NA, const arma::mat& values,
                  const double& lambda, const std::string& type,
                  const double& svd_tol, const double& loss_tuning,
                  const double& delta, double mu, const double& conv_tol,
                  const arma::uword& max_iter,
                  // output to be returned through arguments
                  arma::mat& L, arma::mat& Z, arma::mat& Theta,
                  double& objective, bool& converged,
                  arma::uword& nb_iter) {

  // initializations
  objective = R_PosInf;
  converged = false;
  nb_iter = 0;

  // iterate update steps of Z, L, and Theta
  arma::uword rank, i, j, k, replace_i, nb_values = values.n_rows, which_min;
  arma::mat U, V, L_minus_Z;
  arma::vec d;
  double nuclear_norm, tmp, objective_step2, objective_step2_min,
  loss_norm, loss, loss_min, previous_objective, change;
  while (!converged && nb_iter < max_iter) {

    // step 1: update Z keeping L fixed
    // soft-thresholding of L + 1/mu * Theta

    // compute SVD
    // TODO: add option to use SVD-ALS or softImpute-ALS instead of the SVD
    arma::svd(U, d, V, L + Theta/mu);
    // adjust singular values for our parametrization
    d -= lambda / mu;
    // compute rank of soft-thresholded singular values
    nuclear_norm = 0;
    rank = 0;
    for (j = 0; j < d.n_elem; j++) {
      if (d(j) > svd_tol) {
        nuclear_norm += d(j);
        rank++;
      }
    }
    // initialize Z with zeros
    Z.zeros();
    // if we have a nonzero soft-thresholded singular value, update Z from
    // the soft-thresholded SVD
    if (rank > 0) {
      // efficient implementation of matrix multiplications without copying parts
      for (i = 0; i < n; i++) {
        for (j = 0; j < p; j++) {
          // add contributions from the different rows of U and columns of V'
          for (k = 0; k < rank; k++) Z(i, j) += U(i, k) * d(k) * V(j, k);
        }
      }
    }

    // step 2: update L keeping Z fixed
    // separable problem for missing and observed values in X

    // update L for cells with missing values in X
    for (i = 0; i < ind_NA.n_elem; i++) {
      // index of current cell to be updated
      replace_i = ind_NA(i);
      // column of current cell to be updated
      // TODO: this should be avoided by passing down indices (i, j) 
      // for which cells are missing
      j = replace_i / n;
      // save some computation time in loop
      tmp = -Z(replace_i) + Theta(replace_i)/mu;
      // initialize the minimum
      which_min = 0;
      objective_step2_min = R_PosInf;
      // loop over the different values and choose the one that minimizes the
      // objective function
      for (k = 0; k < nb_values; k++) {
        // The paper says to take the argmin of the squared expression. But
        // this is equivalent to taking the argmin of the absolute value,
        // which is faster to compute.
        objective_step2 = std::abs(values(k, j) + tmp);
        if (objective_step2 < objective_step2_min) {
          which_min = k;
          objective_step2_min = objective_step2;
        }
      }
      // update the element of L with the argmin of the objective function
      L(replace_i) = values(which_min, j);
    }

    // update L for cells with observed values in X
    loss_norm = 0;
    for (i = 0; i < ind_not_NA.n_elem; i++) {
      // index of current cell to be updated
      replace_i = ind_not_NA(i);
      // column of current cell to be updated
      // TODO: this should be avoided by passing down indices (i, j) 
      // for which cells are missing
      j = replace_i / n;
      // save some computation time in loop
      tmp = -Z(replace_i) + Theta(replace_i)/mu;
      // initialize the minimum
      which_min = 0;
      objective_step2_min = R_PosInf;
      // loop over the different values and choose the one that minimizes the
      // objective function
      for (k = 0; k < nb_values; k++) {
        loss = bounded(values(k, j) - X(replace_i), loss_tuning);
        objective_step2 = loss + mu * std::pow(values(k, j) + tmp, 2.0)/2.0;
        if (objective_step2 < objective_step2_min) {
          which_min = k;
          loss_min = loss;
          objective_step2_min = objective_step2;
        }
      }
      // update the element of L with the argmin of the objective function
      L(replace_i) = values(which_min, j);
      // update the norm given by loss function
      loss_norm += loss_min;
    }

    // step 3: update Lagrange multiplier Theta and parameter mu
    L_minus_Z = L - Z;
    Theta = Theta + mu * L_minus_Z;
    mu = delta * mu;

    // update iteration counter
    nb_iter++;
    // update objective function for convergence criterion
    previous_objective = objective;
    objective = loss_norm + lambda * nuclear_norm + arma::dot(Theta, L_minus_Z) +
      mu * arma::norm(L_minus_Z, "fro") / 2.0;
    // compute relative change and check convergence
    if (nb_iter > 1) {
      // we can't compute relative change in the first iteration since the
      // objective function is initialized with infinity
      change = std::abs((objective - previous_objective) / previous_objective);
      converged = change < conv_tol;
    }

  }

}


// ----------------------------------------------------
// function to be called from R without starting values
// ----------------------------------------------------

// [[Rcpp::export]]
Rcpp::List rdmc_cpp(const arma::mat& X, const arma::umat& is_NA,
                    const arma::mat& values, const arma::vec& lambda,
                    const std::string& type, const std::string& loss, 
                    const double& svd_tol, const double& loss_tuning, 
                    const double& delta, double mu, const double& conv_tol, 
                    const arma::uword& max_iter) {
  
  // extract number of rows and columns
  arma::uword n = X.n_rows, p = X.n_cols;
  
  // initialize L and Z
  arma::mat L = X, Z = X;
  // initialize missing values with columnwise median
  arma::vec x_j;
  arma::uvec is_NA_j, keep, replace;
  double median_j;
  arma::uword i, j, replace_i;
  for (j = 0; j < p; j++) {
    x_j = X.unsafe_col(j);
    is_NA_j = is_NA.unsafe_col(j);
    keep = arma::find(is_NA_j == 0);
    median_j = arma::median(x_j.elem(keep));
    replace = arma::find(is_NA_j == 1);
    for (i = 0; i < replace.n_elem; i++) {
      replace_i = replace(i);
      L(replace_i, j) = median_j;
      Z(replace_i, j) = median_j;
    }
  }
  
  // initialize Theta
  arma::mat Theta(n, p);
  Theta.zeros();
  
  // find indices of missing and nonmissing cells of X
  arma::uvec ind_NA = find(is_NA == 1);
  arma::uvec ind_not_NA = find(is_NA == 0);
  
  // initialize variables related to convergence
  // (to be updated by workhorse function)
  double objective;
  bool converged;
  arma::uword nb_iter;
  
  // different behavior depending on whether we have one value of the
  // regularization parameter lambda or multiple values
  if (lambda.n_elem == 1) {
    
    // call workhorse function with initial values
    if (loss == "pseudo_huber") {
      rdmc_pseudo_huber(X, n, p, ind_NA, ind_not_NA, values, lambda(0), type,
                        svd_tol, loss_tuning, delta, mu, conv_tol, max_iter, 
                        L, Z, Theta, objective, converged, nb_iter);
    } else if (loss == "absolute") {
      rdmc_absolute(X, n, p, ind_NA, ind_not_NA, values, lambda(0), type,
                    svd_tol, delta, mu, conv_tol, max_iter, L, Z, Theta,
                    objective, converged, nb_iter);
    } else {
      rdmc_bounded(X, n, p, ind_NA, ind_not_NA, values, lambda(0), type,
                   svd_tol, loss_tuning, delta, mu, conv_tol, max_iter,
                   L, Z, Theta, objective, converged, nb_iter);
    }
    // return list of results
    return Rcpp::List::create(Rcpp::Named("lambda") = lambda(0),
                              Rcpp::Named("L") = L,
                              Rcpp::Named("Z") = Z,
                              Rcpp::Named("Theta") = Theta,
                              Rcpp::Named("objective") = objective,
                              Rcpp::Named("converged") = converged,
                              Rcpp::Named("nb_iter") = nb_iter);
    
  } else {
    
    // loop over values of the regularization parameter lambda
    arma::uword l;
    Rcpp::List out(lambda.n_elem);
    for (l = 0; l < lambda.n_elem; l++) {
      // call workhorse function with starting values: note that solutions
      // for previous value of lambda are used as starting values
      if (loss == "pseudo_huber") {
        rdmc_pseudo_huber(X, n, p, ind_NA, ind_not_NA, values, lambda(l), type,
                          svd_tol, loss_tuning, delta, mu, conv_tol, max_iter, 
                          L, Z, Theta, objective, converged, nb_iter);
      } else if (loss == "absolute") {
        rdmc_absolute(X, n, p, ind_NA, ind_not_NA, values, lambda(l), type,
                      svd_tol, delta, mu, conv_tol, max_iter, L, Z, Theta,
                      objective, converged, nb_iter);
      } else {
        rdmc_bounded(X, n, p, ind_NA, ind_not_NA, values, lambda(l), type,
                     svd_tol, loss_tuning, delta, mu, conv_tol, max_iter,
                     L, Z, Theta, objective, converged, nb_iter);
      }
      // add list of results for current value of lambda: note that a copy of
      // the objects that are stored in the list so that they are not modified
      // in future calls to rdmc_cpp()
      out[l] = Rcpp::List::create(Rcpp::Named("lambda") = lambda(l),
                                  Rcpp::Named("L") = L,
                                  Rcpp::Named("Z") = Z,
                                  Rcpp::Named("Theta") = Theta,
                                  Rcpp::Named("objective") = objective,
                                  Rcpp::Named("converged") = converged,
                                  Rcpp::Named("nb_iter") = nb_iter);
    }
    
    // return list of results for all values of lambda
    return out;
    
  }
  
}


// [[Rcpp::export]]
Rcpp::List rdmc_rank_max_cpp(const arma::mat& X, const arma::umat& is_NA,
                             const arma::mat& values, const arma::vec& lambda, 
                             const arma::uword& rank_max, const std::string& type, 
                             const std::string& loss, const double& svd_tol, 
                             const double& loss_tuning, const double& delta, 
                             double mu, const double& conv_tol, 
                             const arma::uword& max_iter) {
  
  // extract number of rows and columns
  arma::uword n = X.n_rows, p = X.n_cols;
  
  // initialize L and Z
  arma::mat L = X, Z = X;
  // initialize missing values with columnwise median
  arma::vec x_j;
  arma::uvec is_NA_j, keep, replace;
  double median_j;
  arma::uword i, j, replace_i;
  for (j = 0; j < p; j++) {
    x_j = X.unsafe_col(j);
    is_NA_j = is_NA.unsafe_col(j);
    keep = arma::find(is_NA_j == 0);
    median_j = arma::median(x_j.elem(keep));
    replace = arma::find(is_NA_j == 1);
    for (i = 0; i < replace.n_elem; i++) {
      replace_i = replace(i);
      L(replace_i, j) = median_j;
      Z(replace_i, j) = median_j;
    }
  }
  
  // initialize Theta
  arma::mat Theta(n, p);
  Theta.zeros();
  
  // find indices of missing and nonmissing cells of X
  arma::uvec ind_NA = find(is_NA == 1);
  arma::uvec ind_not_NA = find(is_NA == 0);
  
  // initialize variables related to convergence
  // (to be updated by workhorse function)
  double objective;
  bool converged;
  arma::uword nb_iter;
  
  // different behavior depending on whether we have one value of the
  // regularization parameter lambda or multiple values
  if (lambda.n_elem == 1) {
    
    // call workhorse function with initial values
    if (loss == "pseudo_huber") {
      rdmc_pseudo_huber_rank_max(X, n, p, ind_NA, ind_not_NA, values, 
                                 lambda(0), rank_max, type, svd_tol, 
                                 loss_tuning, delta, mu, conv_tol, max_iter, 
                                 L, Z, Theta, objective, converged, nb_iter);
    // } else if (loss == "absolute") {
    //   rdmc_absolute(X, n, p, ind_NA, ind_not_NA, values, lambda(0), type,
    //                 svd_tol, delta, mu, conv_tol, max_iter, L, Z, Theta,
    //                 objective, converged, nb_iter);
    // } else {
    //   rdmc_bounded(X, n, p, ind_NA, ind_not_NA, values, lambda(0), type,
    //                svd_tol, loss_tuning, delta, mu, conv_tol, max_iter,
    //                L, Z, Theta, objective, converged, nb_iter);
    }
    // return list of results
    return Rcpp::List::create(Rcpp::Named("lambda") = lambda(0),
                              Rcpp::Named("L") = L,
                              Rcpp::Named("Z") = Z,
                              Rcpp::Named("Theta") = Theta,
                              Rcpp::Named("objective") = objective,
                              Rcpp::Named("converged") = converged,
                              Rcpp::Named("nb_iter") = nb_iter);
    
  } else {
    
    // loop over values of the regularization parameter lambda
    arma::uword l;
    Rcpp::List out(lambda.n_elem);
    for (l = 0; l < lambda.n_elem; l++) {
      // call workhorse function with starting values: note that solutions
      // for previous value of lambda are used as starting values
      if (loss == "pseudo_huber") {
        rdmc_pseudo_huber_rank_max(X, n, p, ind_NA, ind_not_NA, values, 
                                   lambda(l), rank_max, type, svd_tol, 
                                   loss_tuning, delta, mu, conv_tol, max_iter, 
                                   L, Z, Theta, objective, converged, nb_iter);
      // } else if (loss == "absolute") {
      //   rdmc_absolute(X, n, p, ind_NA, ind_not_NA, values, lambda(l), type,
      //                 svd_tol, delta, mu, conv_tol, max_iter, L, Z, Theta,
      //                 objective, converged, nb_iter);
      // } else {
      //   rdmc_bounded(X, n, p, ind_NA, ind_not_NA, values, lambda(l), type,
      //                svd_tol, loss_tuning, delta, mu, conv_tol, max_iter,
      //                L, Z, Theta, objective, converged, nb_iter);
      }
      // add list of results for current value of lambda: note that a copy of
      // the objects that are stored in the list so that they are not modified
      // in future calls to rdmc_cpp()
      out[l] = Rcpp::List::create(Rcpp::Named("lambda") = lambda(l),
                                  Rcpp::Named("L") = L,
                                  Rcpp::Named("Z") = Z,
                                  Rcpp::Named("Theta") = Theta,
                                  Rcpp::Named("objective") = objective,
                                  Rcpp::Named("converged") = converged,
                                  Rcpp::Named("nb_iter") = nb_iter);
    }
    
    // return list of results for all values of lambda
    return out;
    
  }
  
}


// ------------------------------------------------------------------------
// function to be called from R with starting values
// (only called after tuning the penalty parameter with the optimal lambda)
// ------------------------------------------------------------------------

// Andreas: I had the following comment as a separate line before the 
// corresponding arguments, but that causes an error during compilation!
// starting values (arguments L, Z, Theta) are not passed as a reference 
// so that a copy is created that can safely be modified
// [[Rcpp::export]]
Rcpp::List rdmc_opt_cpp(const arma::mat& X, const arma::umat& is_NA,
                        const arma::mat& values, const double& lambda,
                        const std::string& type, const std::string& loss, 
                        const double& svd_tol, const double& loss_tuning, 
                        const double& delta, double mu, const double& conv_tol, 
                        const arma::uword& max_iter, arma::mat L, arma::mat Z, 
                        arma::mat Theta) {
  
  // extract number of rows and columns
  arma::uword n = X.n_rows, p = X.n_cols;
  
  // find indices of missing and nonmissing cells of X
  arma::uvec ind_NA = find(is_NA == 1);
  arma::uvec ind_not_NA = find(is_NA == 0);
  
  // initialize variables related to convergence
  // (to be updated by workhorse function)
  double objective;
  bool converged;
  arma::uword nb_iter;
  
  // call workhorse function with initial values
  if (loss == "pseudo_huber") {
    rdmc_pseudo_huber(X, n, p, ind_NA, ind_not_NA, values, lambda, type,
                      svd_tol, loss_tuning, delta, mu, conv_tol, max_iter,
                      L, Z, Theta, objective, converged, nb_iter);
  } else if (loss == "absolute") {
    rdmc_absolute(X, n, p, ind_NA, ind_not_NA, values, lambda, type,
                  svd_tol, delta, mu, conv_tol, max_iter, L, Z, Theta,
                  objective, converged, nb_iter);
  } else {
    rdmc_bounded(X, n, p, ind_NA, ind_not_NA, values, lambda, type,
                 svd_tol, loss_tuning, delta, mu, conv_tol, max_iter,
                 L, Z, Theta, objective, converged, nb_iter);
  }
  
  // return list of results
  return Rcpp::List::create(Rcpp::Named("lambda") = lambda,
                            Rcpp::Named("L") = L,
                            Rcpp::Named("Z") = Z,
                            Rcpp::Named("Theta") = Theta,
                            Rcpp::Named("objective") = objective,
                            Rcpp::Named("converged") = converged,
                            Rcpp::Named("nb_iter") = nb_iter);
  
}





// [[Rcpp::export]]
Rcpp::List rdmc_opt_rank_max_cpp(const arma::mat& X, const arma::umat& is_NA,
                                 const arma::mat& values, const double& lambda,
                                 const arma::uword& rank_max, const std::string& type, 
                                 const std::string& loss, const double& svd_tol, 
                                 const double& loss_tuning, const double& delta, 
                                 double mu, const double& conv_tol, 
                                 const arma::uword& max_iter, arma::mat L, 
                                 arma::mat Z, arma::mat Theta) {
  
  // extract number of rows and columns
  arma::uword n = X.n_rows, p = X.n_cols;
  
  // find indices of missing and nonmissing cells of X
  arma::uvec ind_NA = find(is_NA == 1);
  arma::uvec ind_not_NA = find(is_NA == 0);
  
  // initialize variables related to convergence
  // (to be updated by workhorse function)
  double objective;
  bool converged;
  arma::uword nb_iter;
  
  // call workhorse function with initial values
  if (loss == "pseudo_huber") {
    rdmc_pseudo_huber_rank_max(X, n, p, ind_NA, ind_not_NA, values, lambda, 
                               rank_max, type, svd_tol, loss_tuning, delta, mu, 
                               conv_tol, max_iter, L, Z, Theta, objective, 
                               converged, nb_iter);
  // } else if (loss == "absolute") {
  //   rdmc_absolute(X, n, p, ind_NA, ind_not_NA, values, lambda, type,
  //                 svd_tol, delta, mu, conv_tol, max_iter, L, Z, Theta,
  //                 objective, converged, nb_iter);
  // } else {
  //   rdmc_bounded(X, n, p, ind_NA, ind_not_NA, values, lambda, type,
  //                svd_tol, loss_tuning, delta, mu, conv_tol, max_iter,
  //                L, Z, Theta, objective, converged, nb_iter);
  }
  
  // return list of results
  return Rcpp::List::create(Rcpp::Named("lambda") = lambda,
                            Rcpp::Named("L") = L,
                            Rcpp::Named("Z") = Z,
                            Rcpp::Named("Theta") = Theta,
                            Rcpp::Named("objective") = objective,
                            Rcpp::Named("converged") = converged,
                            Rcpp::Named("nb_iter") = nb_iter);
  
}






// -------------------------------------------------------
// previous version that only implements pseudo-Huber loss
// -------------------------------------------------------

// // workhorse function for a single value of the regularization parameter lambda
// void rdmc_cpp(const arma::mat& X, const arma::uword& n, const arma::uword& p,
//               const arma::uvec& ind_NA, const arma::uvec& ind_not_NA,
//               const arma::vec& values, const double& lambda,
//               const char& type, const double& svd_tol, const double& delta,
//               double mu, const double& conv_tol, const arma::uword& max_iter,
//               // output to be returned through arguments
//               arma::mat& L, arma::mat& Z, arma::mat& Theta,
//               double& objective, bool& converged, arma::uword& nb_iter) {
// 
//   // initializations
//   objective = R_PosInf;
//   converged = false;
//   nb_iter = 0;
// 
//   // iterate update steps of Z, L, and Theta
//   arma::uword rank, i, j, k, replace_i, nb_values = values.n_elem, which_min;
//   arma::mat U, V, L_minus_Z;
//   arma::vec d;
//   double nuclear_norm, tmp, objective_step2, objective_step2_min,
//     loss_norm, loss, loss_min, previous_objective, change;
//   while (!converged && nb_iter < max_iter) {
// 
//     // step 1: update Z keeping L fixed
//     // soft-thresholding of L + 1/mu * Theta
// 
//     // compute SVD
//     // TODO: add option to use SVD-ALS or softImpute-ALS instead of the SVD
//     arma::svd(U, d, V, L + Theta/mu);
//     // adjust singular values for our parametrization
//     d -= lambda / mu;
//     // compute rank of soft-thresholded singular values
//     nuclear_norm = 0;
//     rank = 0;
//     for (j = 0; j < d.n_elem; j++) {
//       if (d(j) > svd_tol) {
//         nuclear_norm += d(j);
//         rank++;
//       }
//     }
//     // initialize Z with zeros
//     Z.zeros();
//     // if we have a nonzero soft-thresholded singular value, update Z from
//     // the soft-thresholded SVD
//     if (rank > 0) {
//       // efficient implementation of matrix multiplications without copying parts
//       for (i = 0; i < n; i++) {
//         for (j = 0; j < p; j++) {
//           // add contributions from the different rows of U and columns of V'
//           for (k = 0; k < rank; k++) Z(i, j) += U(i, k) * d(k) * V(j, k);
//         }
//       }
//     }
// 
//     // step 2: update L keeping Z fixed
//     // separable problem for missing and observed values in X
// 
//     // update L for cells with missing values in X
//     for (i = 0; i < ind_NA.n_elem; i++) {
//       // index of current cell to be updated
//       replace_i = ind_NA(i);
//       // save some computation time in loop
//       tmp = -Z(replace_i) + Theta(replace_i)/mu;
//       // initialize the minimum
//       which_min = 0;
//       objective_step2_min = R_PosInf;
//       // loop over the different values and choose the one that minimizes the
//       // objective function
//       for (k = 0; k < nb_values; k++) {
//         // The paper says to take the argmin of the squared expression. But
//         // this is equivalent to taking the argmin of the absolute value,
//         // which is faster to compute.
//         objective_step2 = std::abs(values(k) + tmp);
//         if (objective_step2 < objective_step2_min) {
//           which_min = k;
//           objective_step2_min = objective_step2;
//         }
//       }
//       // update the element of L with the argmin of the objective function
//       L(replace_i) = values(which_min);
//     }
// 
//     // update L for cells with observed values in X
//     loss_norm = 0;
//     for (i = 0; i < ind_not_NA.n_elem; i++) {
//       // index of current cell to be updated
//       replace_i = ind_not_NA(i);
//       // save some computation time in loop
//       tmp = -Z(replace_i) + Theta(replace_i)/mu;
//       // initialize the minimum
//       which_min = 0;
//       objective_step2_min = R_PosInf;
//       // loop over the different values and choose the one that minimizes the
//       // objective function
//       for (k = 0; k < nb_values; k++) {
//         // loss = std::abs(values(k) - X(replace_i));
//         loss = pseudo_huber(values(k) - X(replace_i), 1.0);
//         objective_step2 = loss + mu * std::pow(values(k) + tmp, 2.0)/2.0;
//         if (objective_step2 < objective_step2_min) {
//           which_min = k;
//           loss_min = loss;
//           objective_step2_min = objective_step2;
//         }
//       }
//       // update the element of L with the argmin of the objective function
//       L(replace_i) = values(which_min);
//       // update the norm given by loss function
//       loss_norm += loss_min;
//     }
// 
//     // step 3: update Lagrange multiplier Theta and parameter mu
//     L_minus_Z = L - Z;
//     Theta = Theta + mu * L_minus_Z;
//     mu = delta * mu;
// 
//     // update iteration counter
//     nb_iter++;
//     // update objective function for convergence criterion
//     previous_objective = objective;
//     objective = loss_norm + lambda * nuclear_norm + arma::dot(Theta, L_minus_Z) +
//       mu * arma::norm(L_minus_Z, "fro") / 2.0;
//     // compute relative change and check convergence
//     if (nb_iter > 1) {
//       // we can't compute relative change in the first iteration since the
//       // objective function is initialized with infinity
//       change = std::abs((objective - previous_objective) / previous_objective);
//       converged = change < conv_tol;
//     }
// 
//   }
// 
// }
// 
// 
// // function to be called from R without starting values
// // [[Rcpp::export]]
// Rcpp::List rdmc_cpp(const arma::mat& X, const arma::umat& is_NA,
//                     const arma::vec& values, const arma::vec& lambda,
//                     const char& type, const double& svd_tol,
//                     const double& delta, double mu, const double& conv_tol,
//                     const arma::uword& max_iter) {
// 
//   // extract number of rows and columns
//   arma::uword n = X.n_rows, p = X.n_cols;
// 
//   // initialize L and Z
//   arma::mat L = X, Z = X;
//   // initialize missing values with columnwise median
//   arma::vec x_j;
//   arma::uvec is_NA_j, keep, replace;
//   double median_j;
//   arma::uword i, j, replace_i;
//   for (j = 0; j < p; j++) {
//     x_j = X.unsafe_col(j);
//     is_NA_j = is_NA.unsafe_col(j);
//     keep = arma::find(is_NA_j == 0);
//     median_j = arma::median(x_j.elem(keep));
//     replace = arma::find(is_NA_j == 1);
//     for (i = 0; i < replace.n_elem; i++) {
//       replace_i = replace(i);
//       L(replace_i, j) = median_j;
//       Z(replace_i, j) = median_j;
//     }
//   }
// 
//   // initialize Theta
//   arma::mat Theta(n, p);
//   Theta.zeros();
// 
//   // find indices of missing and nonmissing cells of X
//   arma::uvec ind_NA = find(is_NA == 1);
//   arma::uvec ind_not_NA = find(is_NA == 0);
// 
//   // initialize variables related to convergence
//   // (to be updated by workhorse function)
//   double objective;
//   bool converged;
//   arma::uword nb_iter;
// 
//   // different behavior depending on whether we have one value of the
//   // regularization parameter lambda or multiple values
//   if (lambda.n_elem == 1) {
// 
//     // call workhorse function with initial values
//     rdmc_cpp(X, n, p, ind_NA, ind_not_NA, values, lambda(0), type,
//              svd_tol, delta, mu, conv_tol, max_iter, L, Z, Theta,
//              objective, converged, nb_iter);
//     // return list of results
//     return Rcpp::List::create(Rcpp::Named("lambda") = lambda(0),
//                               Rcpp::Named("L") = L,
//                               Rcpp::Named("Z") = Z,
//                               Rcpp::Named("Theta") = Theta,
//                               Rcpp::Named("objective") = objective,
//                               Rcpp::Named("converged") = converged,
//                               Rcpp::Named("nb_iter") = nb_iter);
// 
//   } else {
// 
//     // loop over values of the regularization parameter lambda
//     arma::uword l;
//     Rcpp::List out(lambda.n_elem);
//     for (l = 0; l < lambda.n_elem; l++) {
//       // call workhorse function with starting values: note that solutions
//       // for previous value of lambda are used as starting values
//       rdmc_cpp(X, n, p, ind_NA, ind_not_NA, values, lambda(l), type,
//                svd_tol, delta, mu, conv_tol, max_iter, L, Z, Theta,
//                objective, converged, nb_iter);
//       // add list of results for current value of lambda: note that a copy of
//       // the objects that are stored in the list so that they are not modified
//       // in future calls to rdmc_cpp()
//       out[l] = Rcpp::List::create(Rcpp::Named("lambda") = lambda(l),
//                                   Rcpp::Named("L") = L,
//                                   Rcpp::Named("Z") = Z,
//                                   Rcpp::Named("Theta") = Theta,
//                                   Rcpp::Named("objective") = objective,
//                                   Rcpp::Named("converged") = converged,
//                                   Rcpp::Named("nb_iter") = nb_iter);
//     }
// 
//     // return list of results for all values of lambda
//     return out;
// 
//   }
// 
// }


// // function to be called from R with starting values
// // [[Rcpp::export]]
// Rcpp::List rdmc_cpp(const arma::mat& X, const arma::umat& is_NA,
//                     const arma::uword& nb_cat, const double& lambda,
//                     const char& type, const double& svd_tol,
//                     const double& delta, double mu, const double& conv_tol,
//                     const arma::uword& max_iter, 
//                     // starting values are not passed as a reference so that
//                     // a copy is created that can safely be modified
//                     const arma::mat L, const arma::mat Z, const arma::mat Theta) {
// 
//   // find indices of missing and nonmissing cells of X
//   arma::uvec ind_NA = find(is_NA == 1);
//   arma::uvec ind_not_NA = find(is_NA == 0);
// 
//   // call workhorse function with starting values
//   double objective;
//   bool converged;
//   arma::uword nb_iter;
//   rdmc_cpp(X, n, p, ind_NA, ind_not_NA, nb_cat, lambda, type, svd_tol, delta,
//            mu, conv_tol, max_iter, L, Z, Theta, objective, converged, nb_iter);
// 
//   // return list of results
//   return Rcpp::List::create(Rcpp::Named("lambda") = lambda(l),
//                             Rcpp::Named("L") = L,
//                             Rcpp::Named("Z") = Z,
//                             Rcpp::Named("Theta") = Theta,
//                             Rcpp::Named("objective") = objective,
//                             Rcpp::Named("converged") = converged,
//                             Rcpp::Named("nb_iter") = nb_iter);
// 
// }
