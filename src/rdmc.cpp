/*
 * Author: Andreas Alfons
 *         Erasmus Universiteit Rotterdam
 * 
 * based on R code by Aurore Archimbaud
 */

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>


// Pseudo-Huber loss
double pseudo_huber(const double& x, const double& delta) {
  return std::pow(delta, 2.0) * (sqrt(1 + std::pow(x/delta, 2.0)) - 1);
}


// workhorse function for a single value of the regularization parameter lambda
void rdmc_cpp(const arma::mat& X, const arma::uword& n, const arma::uword& p, 
              const arma::uvec& ind_NA, const arma::uvec& ind_not_NA, 
              const arma::vec& values, const double& lambda, 
              const char& type, const double& svd_tol, const double& delta, 
              double mu, const double& conv_tol, const arma::uword& max_iter, 
              // output to be returned through arguments
              arma::mat& L, arma::mat& Z, arma::mat& Theta, 
              double& frobenius, bool& converged, arma::uword& nb_iter) {
  
  // initializations
  frobenius = R_PosInf;
  converged = false;
  nb_iter = 0;
  
  // iterate update steps of Z and L
  arma::uword rank, i, j, k, replace_i, nb_values = values.n_elem;
  arma::mat U, V;
  arma::vec d;
  double tmp, which_min, min, objective;
  while (!converged && nb_iter < max_iter) {
    
    // step 1: update Z keeping L fixed
    // soft-thresholding of L + 1/mu * Theta
    
    // compute SVD
    // TODO: add option to use softImpute-ALS or SVD-ALS instead of the SVD
    arma::svd(U, d, V, L + Theta/mu);
    // adjust singular values for our parametrization
    d -= lambda / mu;
    // compute rank of soft-thresholded singular values
    rank = 0;
    for (j = 0; j < d.n_elem; j++) {
      if (d(j) > svd_tol) rank++;
    }
    if (rank == 0) break;
    // update Z from the soft-thresholded SVD
    // (efficient implementation of matrix multiplications without copying parts)
    for (i = 0; i < n; i++) {
      for (j = 0; j < p; j++) {
        // initialize current matrix element
        Z(i, j) = 0.0;
        // add contributions from the different rows of U and columns of V'
        for (k = 0; k < rank; k++) Z(i, j) += U(i, k) * d(k) * V(j, k);
      }
    }
    
    // step 2: update L keeping Z fixed
    // separable problem for missing and observed values in X
    
    // update L for cells with missing values in X
    for (i = 0; i < ind_NA.n_elem; i++) {
      // index of current cell to be updated
      replace_i = ind_NA(i);
      // save some computation time in loop
      tmp = -Z(replace_i) + Theta(replace_i)/mu;
      // initialize the minimum
      which_min = 0;
      min = R_PosInf;
      // loop over the different values and choose the one that minimizes the 
      // objective function
      for (k = 0; k < nb_values; k++) {
        // The paper says to take the argmin of the squared expression. But 
        // this is equivalent to taking the argmin of the absolute value, 
        // which is faster to compute.
        objective = std::abs(values(k) + tmp);
        if (objective < min) {
          which_min = k;
          min = objective;
        }
      }
      // update the element of L with the argmin of the objective function
      L(replace_i) = values(which_min);
    }
    
    // update L for cells with observed values in X
    for (i = 0; i < ind_not_NA.n_elem; i++) {
      // index of current cell to be updated
      replace_i = ind_not_NA(i);
      // save some computation time in loop
      tmp = -Z(replace_i) + Theta(replace_i)/mu;
      // initialize the minimum
      which_min = 0;
      min = R_PosInf;
      // loop over the different values and choose the one that minimizes the 
      // objective function
      for (k = 0; k < nb_values; k++) {
        // objective = std::abs(values(k) - X(replace_i)) + mu * std::pow(values(k) + tmp, 2.0)/2.0;
        objective = pseudo_huber(values(k) - X(replace_i), 1.0) + mu * std::pow(values(k) + tmp, 2.0)/2.0;
        if (objective < min) {
          which_min = k;
          min = objective;
        }
      }
      // update the element of L with the argmin of the objective function
      L(replace_i) = values(which_min);
    }
    
    // step 3: update Lagrange multiplier Theta and parameter mu
    Theta = Theta + mu * (L - Z);
    mu = delta * mu;
    
    // update Frobenius norm (for convergence criterion) and iteration counter
    frobenius = arma::norm(L - Z, "fro");
    converged = frobenius < conv_tol;
    nb_iter++;
    
  }
  
}


// function to be called from R without starting values
// [[Rcpp::export]]
Rcpp::List rdmc_cpp(const arma::mat& X, const arma::umat& is_NA, 
                    const arma::vec& values, const arma::vec& lambda, 
                    const char& type, const double& svd_tol, 
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
  double frobenius;
  bool converged;
  arma::uword nb_iter;
  
  // different behavior depending on whether we have one value of the 
  // regularization parameter lambda or mulitple values
  if (lambda.n_elem == 1) {
    
    // reset Theta to zeros
    // Theta.zeros();
    // call workhorse function with initial values
    rdmc_cpp(X, n, p, ind_NA, ind_not_NA, values, lambda(0), type, 
             svd_tol, delta, mu, conv_tol, max_iter, L, Z, Theta, 
             frobenius, converged, nb_iter);
    // return list of results
    return Rcpp::List::create(Rcpp::Named("lambda") = lambda(0), 
                              Rcpp::Named("L") = L, 
                              Rcpp::Named("Z") = Z, 
                              Rcpp::Named("Theta") = Theta, 
                              Rcpp::Named("frobenius") = frobenius,
                              Rcpp::Named("converged") = converged, 
                              Rcpp::Named("nb_iter") = nb_iter);
    
  } else {
    
    // loop over values of the regularization parameter lambda
    arma::uword l;
    Rcpp::List out(lambda.n_elem);
    for (l = 0; l < lambda.n_elem; l++) {
      // reset Theta to zeros
      // Theta.zeros();
      // call workhorse function with starting values: note that solutions
      // for previous value of lambda are used as starting values
      rdmc_cpp(X, n, p, ind_NA, ind_not_NA, values, lambda(l), type,
               svd_tol, delta, mu, conv_tol, max_iter, L, Z, Theta,
               frobenius, converged, nb_iter);
      // add list of results for current value of lambda: note that a copy of
      // the objects that are stored in the list so that they are not modified
      // in future calls to rdmc_cpp()
      out[l] = Rcpp::List::create(Rcpp::Named("lambda") = lambda(l),
                                  Rcpp::Named("L") = L,
                                  Rcpp::Named("Z") = Z,
                                  Rcpp::Named("Theta") = Theta,
                                  Rcpp::Named("frobenius") = frobenius,
                                  Rcpp::Named("converged") = converged,
                                  Rcpp::Named("nb_iter") = nb_iter);
    }
    
    // return list of results for all values of lambda
    return out;
    
  }
  
}


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
//   double frobenius;
//   bool converged;
//   arma::uword nb_iter;
//   rdmc_cpp(X, n, p, ind_NA, ind_not_NA, nb_cat, lambda, type, svd_tol, delta,
//            mu, conv_tol, max_iter, L, Z, Theta, frobenius, converged, nb_iter);
// 
//   // return list of results
//   return Rcpp::List::create(Rcpp::Named("L") = L,
//                             Rcpp::Named("Z") = Z,
//                             Rcpp::Named("Theta") = Theta,
//                             Rcpp::Named("frobenius") = frobenius,
//                             Rcpp::Named("converged") = converged,
//                             Rcpp::Named("nb_iter") = nb_iter);
// 
// }
