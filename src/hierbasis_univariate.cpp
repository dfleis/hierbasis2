// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "helpers.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::sp_mat prox(arma::vec u, arma::mat w) {
  // Implementation of the proximal method outlined in
  // Jenatton et al. (2010) (via Haris et al. 2016).
  //
  // Solves the minimization problem
  //   (1/2) * || u - \beta ||_2^2 + P(\beta),
  // where
  //   P(\beta) = \sum_{i=1}^{p} w[i, j] * l2norm(beta[i:p]),
  // for a column vector  w[,j].
  //
  // Args:
  //    * u: A p-vector of transformed responses (u = U^T %*% y / n).
  //    * w: A matrix of weights used in the penalty function such that
  //         each column corresponds to a fixed penalty \lambda.
  // Returns:
  //    * beta: A p-by-ncol(w) real-valued matrix of solutions
  //            to the proximal gradient descent problem.

  // extract dimensions
  int p    = u.size(); // size of u = number of basis functions
  int nlam = w.n_cols; // number of columns = number of penalties

  // initialize the sparse matrix to store our optima
  arma::sp_mat beta(p, nlam);
  // initialize temporary weights to be used in each iteration
  arma::vec w_tmp;

  for (int j = 0; j < nlam; ++j) {
    // work on j^th column/penalty weights
    w_tmp = w.col(j);
    // temporary beta vector for each iteration of the proximal loop
    arma::vec beta_tmp(u.begin(), u.size(), true);

    // main loop for solving the proximal problem
    for (int i = p - 1; i >= 0; --i) {
      // update estimates
      beta_tmp.subvec(i, p - 1) =
        max(1 - w_tmp(i) / norm(beta_tmp.subvec(i, p - 1)), 0) *
        beta_tmp.subvec(i, p - 1);
    }

    // place the j^th solution in the sparse matrix
    beta.col(j) = beta_tmp;
  }

  return beta;
}

// [[Rcpp::export]]
List solvehierbasis(arma::mat design_mat,
                    arma::vec y,
                    arma::vec ak,
                    arma::mat weights,
                    int n,
                    double lam_min_ratio,
                    int nlam,
                    double max_lambda) {
  // This function returns the argmin of the following function:
  // (0.5) * ||y - X %*% beta||_2^2 + P(beta), where
  // P(\beta) = \sum_{j=1}^{J_n} weights[j, lam_indx] * norm(beta[j:J_n])
  // where j indexes the basis function and lam_indx indexes the different
  // lambda values.
  //
  // Args:
  //    design_mat: A design matrix of size n * J_n, where J_n is number of
  //                basis functions. This matrix should be centered.
  //    y: The centered response vector.
  //    ak: A J_n-vector of weights where ak[j] = j^m - (j-1)^m for a smoothness
  //        level m.
  //    weights: A J_n * nlam matrix which is simply the concatenation [ak,ak,...,ak].
  //    n: The number of observations, equal to length(y).
  //    lam_min_ratio: Ratio of min_lambda to max_lambda.
  //    nlam: Number of lambda values.
  //    max_lambda: A double specifying the maximum lambda, if NA then function
  //                selects a max_lambda.
  // Returns:
  //    beta_hat2: The solution matrix, a J_n * nlam  matrix.
  //    lambdas: Sequence of lambda values used for fitting.

  // Generate the x_mat, this is our orthonormal design matrix.
  arma::mat x_mat, R_mat;
  arma::qr_econ(x_mat, R_mat, design_mat);

  x_mat = x_mat * sqrt(n);
  R_mat  = R_mat / sqrt(n);
  // arma::sp_mat R_mat2 = sp_mat(R_mat / sqrt(n));


  // Note that for a orthogonal design with X^T %*% X = nI the
  // two problems are equivalent:
  // 1. (1/2n)*|Y - X %*% beta|^2 + lambda * Pen(beta),
  // 2. (1/2)*|t(X) %*% Y/n - beta|^2 + lambda * Pen(beta).
  //
  // Thus all we need, is to solve the proximal problem.
  // Define v_temp = t(X) %*% Y/n for the prox problem.

  arma::vec v_temp = x_mat.t() * (y /n);

  // We evaluate a max_lambda value if not specified by user.
  if(R_IsNA(max_lambda)) {
    // Followed by the maximum lambda value.
    max_lambda = max(abs(v_temp) / ak);
  }

  // Generate the full lambda sequence.
  arma::vec lambdas = linspace<vec>(log10(max_lambda),
                                    log10(max_lambda * lam_min_ratio),
                                    nlam);
  lambdas = exp10(lambdas);

  // Generate matrix of weights.
  weights.each_row() %= lambdas.t();

  // Obtain beta.hat, solve prox problem.
  arma::sp_mat beta_hat =  prox(v_temp, weights);

  // Take the inverse of the R_mat to obtain solution on the original scale
  arma::mat beta_hat2 = solve(trimatu(R_mat), mat(beta_hat));

  return List::create(Named("beta") = beta_hat2,
                      Named("lambdas") = lambdas);

}

