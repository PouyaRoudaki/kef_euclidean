#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>  // For parallel computing
#endif
#include "density.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

// Function to compute the marginal log likelihood for a pair of hyperparameters lambda and tau
// [[Rcpp::export]]
double marginal_log_likelihood(
    const arma::mat& centered_kernel_mat_samples,
    const arma::vec& samples,
    const arma::vec& base_measure_weights,
    double dimension,
    const std::string& data_type,
    const arma::vec p_vec,
    double lambda,
    double tau,
    const arma::mat& std_rnorm_matrix,
    int MC_iterations,
    bool parallel_computing = true) {

  int n = centered_kernel_mat_samples.n_rows;

  arma::mat w_sampled(MC_iterations, n); // Make the w_sampled matrix
  //p_vec = arma::vectorise(p_vec);  // Ensure column vector

  // Parallel Monte Carlo weight sampling
#pragma omp parallel for if(parallel_computing) num_threads(omp_get_max_threads())
  for (int i = 0; i < MC_iterations; i++) {
    for (int j = 0; j < n; j++) {
      w_sampled(i, j) = std_rnorm_matrix(i, j) * std::sqrt(p_vec(j) / tau); //std::sqrt(p_vec(j) / tau)
    }
  }

  // Compute densities for sampled weights
  arma::mat densities_for_given_weights(MC_iterations, n);
  bool nan_detected = false;

#pragma omp parallel for if(parallel_computing) num_threads(omp_get_max_threads())
  for (int i = 0; i < MC_iterations; i++) {
    arma::vec dens_row = get_dens_wo_grid(
      centered_kernel_mat_samples, samples, base_measure_weights,
      dimension, data_type, lambda, w_sampled.row(i).t()
    );

    if (!dens_row.is_finite()) {
#pragma omp critical
{
  Rcpp::Rcout << "Error: NaN values in get_dens_wo_grid at iteration " << i << std::endl;
}
      nan_detected = true;
    }

    densities_for_given_weights.row(i) = dens_row.t();
  }

  if (nan_detected) {
    return -arma::datum::inf;
  }

  // Compute the product along each column (equivalent to apply(..., 2, prod) in R)
  arma::vec monte_carlo_likelihood = arma::prod(densities_for_given_weights, 1);


  //Rcpp::Rcout << "p_vec size: " << monte_carlo_likelihood.n_rows << " x " << monte_carlo_likelihood.n_cols << std::endl;

  // Compute the marginal likelihood (mean of all joint conditional densities)
  double marginal_likelihood = arma::mean(monte_carlo_likelihood);

  double log_marginal_likelihood = std::log(marginal_likelihood);

  // Handle -inf case
  if (!std::isfinite(log_marginal_likelihood)) {
    log_marginal_likelihood = -10000;
  }

  //Rcpp::Rcout << "-----------------------------------------------------" << std::endl;
  //Rcpp::Rcout << "lambda: " << lambda << ", tau: " << tau << ", mll: " << log_marginal_likelihood <<std::endl;
  //Rcpp::Rcout << "-----------------------------------------------------" << std::endl;

  // Return the marginal_log_likelihood
  return log_marginal_likelihood;
}
