#ifndef DENSITY_H  // Prevents multiple inclusions
#define DENSITY_H

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// Function declaration (prototype)
arma::vec get_dens_wo_grid(const arma::mat& centered_kernel_mat_samples,
                           const arma::mat& samples,
                           const arma::vec& base_measure_weights,
                           double dimension,
                           const std::string& data_type,
                           double lambda,
                           const arma::vec& weight_vec);

#endif  // DENSITY_H
