#' Compute Base Measure Weights
#'
#' Computes the base measure weights for either 1D (using grid midpoints) or 2D (using Voronoi tessellation).
#'
#' @param sample A numeric vector or matrix of sample points. For 1D, a vector; for 2D, a matrix/data.frame with two columns.
#' @param boundaries A numeric vector (length 2) for 1D defining min/max bounds, or a 2x2 matrix for 2D, each row defining min/max bounds of each axis.
#' @param dimension Integer. Dimensionality of the sample: either 1 or 2.
#'
#' @return A numeric vector of base measure weights.
#' @importFrom deldir deldir
#' @export
#'
#' @examples
#' get_base_measures(sample = runif(100), boundaries = c(0,1), dimension = 1)
#' get_base_measures(sample = matrix(runif(200), ncol=2), boundaries = matrix(c(0,1,0,1), nrow=2), dimension = 2)
get_base_measures <- function(samples, boundaries, data_type = "euclidean",
                              dimension = 1) {
  if (data_type == "euclidean"){

    if (dimension == 1) {
      # Validate input
      if (length(boundaries) != 2) stop("For 1D, 'boundaries' must be a length-2 vector.")

      # Get the sort index and reverse mapping
      sort_index <- order(samples)
      unsort_index <- order(sort_index)  # this gives the positions to map back

      # Sorted samples for computing base measures
      samples_sorted <- samples[sort_index]

      # Compute middle points between sorted samples
      min_x <- boundaries[1]
      max_x <- boundaries[2]
      middle_points <- get_middle_points_grid(min_x, samples_sorted, max_x)

      # Compute interval lengths (base measure on sorted)
      base_measure_weights_sorted <- diff(middle_points)

      # Reorder to match original sample order
      base_measure_weights <- base_measure_weights_sorted[unsort_index]

    } else if (dimension == 2) {
      # Validate input
      if (!is.matrix(boundaries) || dim(boundaries)[1] != 2 || dim(boundaries)[2] != 2)
        stop("For 2D, 'boundaries' must be a 2x2 matrix.")
      if (!is.matrix(samples) || ncol(samples) != 2)
        stop("For 2D, 'samples' must be a matrix with 2 columns.")

      x1 <- samples[, 1]
      x2 <- samples[, 2]
      min_x1 <- boundaries[1, 1]
      max_x1 <- boundaries[1, 2]
      min_x2 <- boundaries[2, 1]
      max_x2 <- boundaries[2, 2]

      # Use Voronoi diagram from deldir
      base_measure_weights <- deldir::deldir(x1, x2,
                                             rw = c(min_x1, max_x1, min_x2, max_x2))$summary$dir.area
    } else {
      stop("Currently, the function only supports 1D or 2D input.")
    }

  } else if(data_type == "order"){
    .Call(`_kefV1_base_measure_order`, samples)
  } else if(data_type == "graph"){
    stop("Graph data detected: this part of the function has not been implemented yet.")
  } else{
    stop("Incompatible data type!")
  }

  return(base_measure_weights)
}


