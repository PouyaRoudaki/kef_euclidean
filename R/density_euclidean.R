#' Compute Unnormalized Density at Sampled Points
#'
#' @param centered_kernel_mat_samples Matrix of centered kernel values for samples.
#' @param lambda Scalar lambda parameter.
#' @param weight_vec Vector of weights.
#' @return Unnormalized density values at sampled points.
#' @export
unnormalised_density_samples <- function(centered_kernel_mat_samples, lambda, weight_vec) {
  .Call(`_kefV2_unnormalised_density_samples`, centered_kernel_mat_samples, lambda, weight_vec)
}

#' Compute Unnormalized Density at Grid Points
#'
#' @param centered_kernel_mat_grids Matrix of centered kernel values between grid and sample points.
#' @param centered_kernel_self_grids Vector of self-kernel values at grid points.
#' @param lambda Scalar lambda parameter.
#' @param weight_vec Vector of weights.
#' @return Unnormalized density values at grid points.
#' @export
unnormalised_density_grids <- function(centered_kernel_mat_grids,
                                       centered_kernel_self_grids, lambda,
                                       weight_vec) {
  .Call(`_kefV2_unnormalised_density_grids`, centered_kernel_mat_grids,
        centered_kernel_self_grids, lambda, weight_vec)
}

#' Compute Normalized Densities at Sampled Points Without a Grid for euclidean data.
#'
#' @param centered_kernel_mat_samples Centered kernel matrix at sampled points.
#' @param samples Sample locations.
#' @param base_measure_weights Vector of base measure weights (used for normalizing constant if dim > 1).
#' @param dimension Dimensionality of the data (1 or higher).
#' @param lambda Scalar lambda parameter.
#' @param weight_vec Vector of weights.
#' @return Normalized density values at sampled points.
#' @export
get_dens_wo_grid_euclidean <- function(centered_kernel_mat_samples, samples,
                                       base_measure_weights, dimension,
                                       #data_type,
                                       lambda, weight_vec) {
  .Call(`_kefV2_get_dens_wo_grid_euclidean`, centered_kernel_mat_samples, samples,
        base_measure_weights, dimension, #data_type,
        lambda, weight_vec)
}

#' Compute Normalized Densities at Sampled and Grid Points for euclidean data.
#'
#' @param centered_kernel_mat_samples Centered kernel matrix at sampled points.
#' @param centered_kernel_mat_grids Centered kernel matrix at grid points.
#' @param centered_kernel_self_grids Self-kernel values at grid points.
#' @param samples Sampled x values.
#' @param grids Grid x values.
#' @param base_measure_weights_grid Vector of base measure weights (used for normalizing constant if dim > 1).
#' @param dimension Dimensionality of the data (1 or higher).
#' @param lambda Scalar lambda parameter.
#' @param weight_vec Vector of weights.
#' @return A list with normalized densities at sampled and grid points.
#' @export
get_dens_euclidean <- function(centered_kernel_mat_samples,
                     centered_kernel_mat_grids,
                     centered_kernel_self_grids,
                     samples,
                     grids,
                     base_measure_weights_grid,
                     dimension,
                     #data_type,
                     lambda,
                     weight_vec) {
  .Call(`_kefV2_get_dens_euclidean`, centered_kernel_mat_samples, centered_kernel_mat_grids,
        centered_kernel_self_grids, samples, grids, base_measure_weights_grid,
        dimension, #data_type,
        lambda, weight_vec)
}




