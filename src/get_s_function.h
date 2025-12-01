#ifndef GET_S_FUNCTION_H
#define GET_S_FUNCTION_H

#include <RcppArmadillo.h>

// Declare the get_s_function function
arma::vec get_s_function(const arma::vec& weight_vec,
                         double lambda,
                         double tau,
                         const arma::mat& centered_kernel_mat_samples,
                         const arma::mat& samples,
                         const arma::vec& base_measure_weights,
                         double dimension,
                         const std::string& data_type,
                         bool prior_var_prob);

#endif // GET_S_FUNCTION_H
