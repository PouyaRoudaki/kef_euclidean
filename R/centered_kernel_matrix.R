#' Compute the Centered Kernel Matrix
#'
#' This function computes the centered fractional Brownian motion (fBM) kernel matrix.
#' @param dimension Integer. Dimensionality of the input space.
#' @param eval_points_1 Matrix or vector. First set of evaluation points (e.g., \eqn{\bx}).
#' @param eval_points_2 Matrix or vector. Second set of evaluation points (e.g., \eqn{\bx'}).
#' @param centering_points Matrix or vector. Grid used for kernel centering (e.g., \eqn{\bz_1, \ldots, \bz_n}).
#' @param hurst_coef Numeric scalar. The Hurst coefficient.
#' @return A numeric matrix of kernel evaluations.
#' @export
centered_kernel_matrix <- function(eval_points_1, eval_points_2, centering_grid,
                                   hurst_coef, dimension, data_type = "euclidean") {
  if(data_type == "euclidean"){
    if (dimension == 1) {
      .Call(`_kefV1_centered_kernel_matrix`, eval_points_1, eval_points_2, centering_grid, hurst_coef)
    } else {
      .Call(`_kefV1_centered_kernel_matrix_hd`, eval_points_1, eval_points_2, centering_grid, hurst_coef)
    }
  }else if(data_type == "order"){
    .Call(`_kefV1_centered_kernel_matrix_order`, eval_points_1, eval_points_2, centering_grid, hurst_coef)
  }else if(data_type == "graph"){
    stop("Graph data detected: this part of the function has not been implemented yet.")
  }else{
    stop("Incompatible data type!")
  }


}
