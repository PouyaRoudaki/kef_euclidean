#' Compute the `s` Function from Weights and Kernel Matrix for Euclidean data
#'
#' Wrapper for the Rcpp function `get_s_function_euclidean` that computes a vector used in density-based optimization.
#'
#' @param weight_vec A numeric vector of weights.
#' @param lambda A scalar regularization parameter.
#' @param tau A scalar penalty parameter.
#' @param centered_kernel_mat_samples Centered kernel matrix at sampled points (n × n).
#' @param samples A numeric matrix of sampled values.
#' @param base_measure_weights A numeric vector of base measure weights (same length as `samples`).
#' @param dimension A scalar indicating the dimension of the data.
#' @param prior_var_prob Logical; if TRUE, the variance of prior is proportional to probability itself and a hyper-parameter, else it is only proportional to a hyper-parameter.
#'
#' @return A numeric vector representing the output of the s-function.
#' @export
get_s_function_euclidean <- function(weight_vec,
                           lambda,
                           tau,
                           centered_kernel_mat_samples,
                           samples,
                           base_measure_weights,
                           dimension,
                           #data_type,
                           prior_var_prob = TRUE) {
  # Force samples to be a numeric matrix (n x d)
  samples_mat <- if (is.matrix(samples)) {
    samples
  } else {
    as.matrix(samples)
  }

  .Call(`_kefV2_get_s_function_euclidean`, weight_vec, lambda, tau,
        centered_kernel_mat_samples, samples,
        base_measure_weights, dimension, #data_type,
        prior_var_prob)
}

