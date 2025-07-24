#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// Compute unnormalized density at sampled points
// [[Rcpp::export]]
arma::vec unnormalised_density_samples(const arma::mat& centered_kernel_mat_samples,
                                       double lambda,
                                       const arma::vec& weight_vec) {
  // Extract diagonal elements for self-kernel correction
  arma::vec diag_vals = centered_kernel_mat_samples.diag();

  // Compute the exponent term for the unnormalized density
  arma::vec exponent_samples = lambda * (centered_kernel_mat_samples.t() * weight_vec - 0.5 * diag_vals);

  // Exponentiate to get unnormalized density
  arma::vec unnorm_dens_vec = arma::exp(exponent_samples);

  return unnorm_dens_vec;
}

// Compute unnormalized density at grid points
// [[Rcpp::export]]
arma::vec unnormalised_density_grids(const arma::mat& centered_kernel_mat_grids,
                                     const arma::vec& centered_kernel_self_grids,
                                     double lambda,
                                     const arma::vec& weight_vec) {
  // Compute the exponent term for the unnormalized density
  arma::vec exponent_grids = lambda * (centered_kernel_mat_grids.t() * weight_vec - 0.5 * centered_kernel_self_grids);

  // Exponentiate to get unnormalized density
  arma::vec unnorm_dens_vec = arma::exp(exponent_grids);

  return unnorm_dens_vec;
}

// Compute normalized densities at sampled points without a grid
// [[Rcpp::export]]
arma::vec get_dens_wo_grid(const arma::mat& centered_kernel_mat_samples,
                           const arma::vec& samples,
                           const arma::vec& base_measure_weights,
                           double dimension,
                           double lambda,
                           const arma::vec& weight_vec) {

  if (arma::any(base_measure_weights < 0)) {
    Rcpp::Rcout << "Warning: Negative base measure weights detected.\n";
  }

  // Compute log-densities (unnormalized log-likelihood contributions)
  arma::vec exponent = lambda * (centered_kernel_mat_samples.t() * weight_vec -
    0.5 * centered_kernel_mat_samples.diag());

  // Handle potential overflow in exponential calculation
  const double EXP_THRESHOLD = 700;
  double max_exponent = arma::max(exponent);
  arma::vec unnorm_density_samples;

  if (max_exponent > EXP_THRESHOLD) {
    // Apply log-stabilization trick: exp(f_x - max_f_x)
    unnorm_density_samples = arma::exp(exponent - max_exponent);
  } else {
    unnorm_density_samples = arma::exp(exponent);
  }

  // Compute normalization constant
  double normalizing_cte;
  if (dimension == 1) {
    // 1D case: use trapezoidal rule for integration
    normalizing_cte = arma::as_scalar(trapz(samples, unnorm_density_samples));
    //normalizing_cte = arma::dot(base_measure_weights, unnorm_density_samples);
  } else {
    // Multivariate case: use dot product with base measure weights
    normalizing_cte = arma::dot(base_measure_weights, unnorm_density_samples);
    //if (normalizing_cte == 0.0){
      //double min_val = arma::min(weight_vec);
      //double max_val = arma::max(weight_vec);
      //double mean_val = arma::mean(weight_vec);
      //int n = weight_vec.size();

      //Rcpp::Rcout << "Summary of weight_vec:\n";
      //Rcpp::Rcout << "Length: " << n << "\n";
      //Rcpp::Rcout << "Min: " << min_val << ", Max: " << max_val << ", Mean: " << mean_val << "\n";
      //normalizing_cte = 1.0;
    //}
  }

  // Check for invalid normalization constant
  if (normalizing_cte == 0.0 || std::isnan(normalizing_cte) || std::isinf(normalizing_cte)) {
    Rcpp::Rcout << "Error: Normalizing constant is zero, NaN, or Inf.\n";
    return arma::vec(samples.n_elem, arma::fill::zeros);  // Return zero vector to avoid NaNs
  }

  arma::vec dens_samples;
  dens_samples = unnorm_density_samples / normalizing_cte;

  if (normalizing_cte < 0) {
    Rcpp::Rcout << "Error: Negative normalizing constant detected. Debug inputs.\n";
  }

  // Return normalized density
  return dens_samples;
}

// Compute normalized densities at sampled and grid points
// [[Rcpp::export]]
Rcpp::List get_dens(const arma::mat& centered_kernel_mat_samples,
                    const arma::mat& centered_kernel_mat_grids,
                    const arma::vec& centered_kernel_self_grids,
                    const arma::vec& samples,
                    const arma::vec& grids,
                    const arma::vec& base_measure_weights,
                    double dimension,
                    double lambda,
                    const arma::vec& weight_vec) {

  // Compute unnormalized density at sampled points
  arma::vec unnorm_density_samples = unnormalised_density_samples(centered_kernel_mat_samples, lambda, weight_vec);

  // Compute unnormalized density at grid points
  arma::vec unnorm_density_grids = unnormalised_density_grids(centered_kernel_mat_grids, centered_kernel_self_grids, lambda, weight_vec);

  // Compute normalization constant
  double normalizing_cte;
  if (dimension == 1) {
    // 1D case: integrate using trapezoidal rule over grid
    //normalizing_cte = arma::as_scalar(trapz(grids, unnorm_density_grids));
    normalizing_cte = arma::dot(base_measure_weights, unnorm_density_samples);
  } else {
    // Higher-dimensional case: weighted sum over base measure weights
    normalizing_cte = arma::dot(base_measure_weights, unnorm_density_samples);
  }

  // Normalize densities
  arma::vec dens_samples_norm = unnorm_density_samples / normalizing_cte;
  arma::vec dens_grid_norm = unnorm_density_grids / normalizing_cte;

  // Return normalized densities
  return Rcpp::List::create(
    Rcpp::Named("samples") = dens_samples_norm,
    Rcpp::Named("grids") = dens_grid_norm,
    Rcpp::Named("unnorm_samples") = unnorm_density_samples,
    Rcpp::Named("unnorm_grids") = unnorm_density_grids,
    Rcpp::Named("norm_cte") = normalizing_cte
  );
}
