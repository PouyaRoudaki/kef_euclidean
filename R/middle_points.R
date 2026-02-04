#' Get Middle Points of Grid (Fast C++ Version)
#'
#' This function calculates the midpoints between adjacent elements of a sample
#' and includes the minimum and maximum grid boundaries.
#'
#' @param min The minimum value of the grid.
#' @param samples A numeric vector representing the sample points within the grid.
#' @param max The maximum value of the grid.
#'
#' @return A numeric vector containing the minimum value, the midpoints of the sample, and the maximum value.
#' @export
#'
#' @examples
#' get_middle_points_grid(0, c(2, 4, 6), 8)
get_middle_points_grid <- function(min, samples, max) {
  .Call(`_kefV2_get_middle_points_grid`, min, samples, max)
}
