#' Compute the `s` Function from Weights and Kernel Matrix
#'
#' Wrapper for the Rcpp function `get_s_function` that computes a vector used in density-based optimization.
#'
#' @param weight_vec A numeric vector of weights.
#' @param lambda A scalar regularization parameter.
#' @param tau A scalar penalty parameter.
#' @param centered_kernel_mat_samples Centered kernel matrix at sampled points (n × n).
#' @param samples A numeric vector of sampled values.
#' @param base_measure_weights A numeric vector of base measure weights (same length as `samples`).
#' @param dimension A scalar indicating the dimensionality of the data.
#'
#' @return A numeric vector representing the output of the s-function.
#' @export
get_s_function <- function(weight_vec,
                           lambda,
                           tau,
                           centered_kernel_mat_samples,
                           samples,
                           base_measure_weights,
                           dimension) {
  .Call(`_kefV1_get_s_function`, weight_vec, lambda, tau,
        centered_kernel_mat_samples, samples,
        base_measure_weights, dimension)
}
