#include <RcppArmadillo.h>
#include "density_euclidean.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// Compute the function s(weight_vec)
// [[Rcpp::export]]
arma::vec get_s_function_euclidean(const arma::vec& weight_vec,
                         double lambda,
                         double tau,
                         const arma::mat& centered_kernel_mat_samples,
                         const arma::mat& samples,
                         const arma::vec& base_measure_weights,
                         double dimension,
                         bool prior_var_prob) {

  // Step 1: Sample size
  double n = centered_kernel_mat_samples.n_rows;
  //Rcout << "[1] Sample size (n): " << n << std::endl;

  // Step 2: Compute densities
  arma::vec dens = get_dens_wo_grid_euclidean(centered_kernel_mat_samples,
                                    samples,
                                    base_measure_weights,
                                    dimension,
                                    lambda,
                                    weight_vec);
  //Rcout << "[2] Densities (dens):\n" << dens.t() << std::endl;

  // Step 3: Density-based probabilities
  arma::vec dens_sample_via_base = dens % base_measure_weights;
  //Rcout << "[3] dens % base_measure_weights:\n" << dens_sample_via_base.t() << std::endl;

  // Step 4: Normalized probability vector
  double total_mass = sum(dens_sample_via_base);
  //Rcout << "[4] Sum of dens_sample_via_base: " << total_mass << std::endl;

  arma::vec prob_sample_via_base = dens_sample_via_base / total_mass;
  //Rcout << "[5] Normalized probabilities (prob_sample_via_base):\n" << prob_sample_via_base.t() << std::endl;

  // Step 5: Column sums of kernel matrix
  arma::rowvec col_sums = sum(centered_kernel_mat_samples, 0);
  //Rcout << "[6] Column sums of kernel matrix (col_sums):\n" << col_sums << std::endl;

  // Step 6: prob_sample_via_base^T * kernel
  arma::rowvec prob_times_kernel = prob_sample_via_base.t() * centered_kernel_mat_samples;
  //Rcout << "[7] prob_sample_via_base^T * kernel (prob_times_kernel):\n" << prob_times_kernel << std::endl;

  // Step 7: Final s computation
  arma::vec s;
  if(prior_var_prob){
    s = lambda * (col_sums.t() - n * prob_times_kernel.t()) -
      tau * (weight_vec / prob_sample_via_base);
  }else{
    s = lambda * (col_sums.t() - n * prob_times_kernel.t()) -
      tau * (weight_vec);
  }

  //Rcout << "[8] Final s vector:\n" << s.t() << std::endl;
  return s;
}


