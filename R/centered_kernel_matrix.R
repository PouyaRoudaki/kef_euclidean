#' Compute a Centered fBM Kernel Matrix
#'
#' Computes the centered fractional Brownian motion (fBM) kernel matrix for
#' Euclidean, order, or graph data.
#'
#' @param eval_points_1 Matrix or numeric vector. First set of evaluation points
#'   (e.g., x).
#' @param eval_points_2 Matrix or numeric vector. Second set of evaluation points
#'   (e.g., x').
#' @param centering_grid Matrix or numeric vector. Grid used for kernel centering
#'   (e.g., z1, ..., zn).
#' @param hurst_coef Numeric scalar. The Hurst coefficient.
#' @param dimension Integer. Dimensionality of the input space when
#'   data_type = "euclidean".
#' @param data_type String specifying the data type. One of "euclidean",
#'   "order", or "graph". The default is "euclidean".
#'
#' @return A numeric matrix of kernel evaluations.
#' @export
centered_kernel_matrix <- function(eval_points_1, eval_points_2, centering_grid,
                                   hurst_coef, dimension, data_type = "euclidean") {
  if (data_type == "euclidean") {
    if (dimension == 1) {
      .Call(`_kefV1_centered_kernel_matrix`, eval_points_1, eval_points_2, centering_grid, hurst_coef)
    } else {
      .Call(`_kefV1_centered_kernel_matrix_hd`, eval_points_1, eval_points_2, centering_grid, hurst_coef)
    }
  } else if (data_type == "order") {
    .Call(`_kefV1_centered_kernel_matrix_order`, eval_points_1, eval_points_2, centering_grid, hurst_coef)
  } else if (data_type == "graph") {
    stop("Graph data detected: this part of the function has not been implemented yet.")
  } else {
    stop("Incompatible data type!")
  }
}
