#' Kernel Density Estimation Using the KEF Method
#'
#' This function estimates density using kernel exponential family (KEF) methods.
#' It constructs a kernel matrix and computes weights using the Barzilai-Borwein
#' optimization method. If the samples size is large, it avoids computational complexity
#' by returning density estimates only at samples points.
#'
#' @param samples A numeric vector representing the observed samples points.
#' @param grids A numeric vector representing the evaluation grids.
#' @param lambda A numeric scalar for the lambda parameter (regularization term).
#' @param tau A numeric scalar for the tau parameter.
#' @param data_type  String specifying the data type. Options are `"euclidean"`, `"order"`, or `"graph"`. The default is `"euclidean"`.
#' @param boundaries A numeric vector (for 1D) or matrix (for multi-dimensional input) specifying the domain boundaries. Optional; if not provided, the domain boundaries are extended by 10% around the samples range.
#' @param prior_var_prob Logical; if TRUE, the variance of prior is proportional to probability itself and a hyper-parameter, else it is only proportional to a hyper-parameter. Default is TRUE.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{weights} - The estimated weight vector.
#'   \item \code{probs_samples} - Density estimates at samples points.
#'   \item \code{probs_grids} - Density estimates at grids points (only if \code{length(samples) < 1000}).
#'   \item \code{time} - The function's execution time.
#' }
#'
#' @details
#' If the samples size is \code{>= 1000}, the function **only returns densities at samples points**
#' to avoid computational complexity. A warning is issued in this case.
#'
#' @export
#'
#' @examples
#' samples <- rnorm(500)
#' grids <- seq(-3, 3, length.out = 100)
#' lambda <- 0.1
#' tau <- 0.5
#' result <- kef(samples = samples, grids = grids,lambda = lambda, tau = tau, data_type="euclidean")
#' plot(grids, result$grids, type = "l", main = "Estimated Density", ylab = "Density")

kef <- function(samples, grids,
                lambda, tau,
                data_type = "euclidean",
                boundaries = NULL,
                prior_var_prob = TRUE) {

  # Start timer
  start_time <- Sys.time()

  if (data_type == "euclidean"){
    message("Euclidean data detected.")
    # Setting default boundaries with 10% padding if the boundaries are not provided.
    if (is.null(boundaries)) {
      if (is.vector(samples)) {
        min_x <- min(samples)
        max_x <- max(samples)
        padding <- 0.1 * (max_x - min_x)
        boundaries <- c(min_x - padding, max_x + padding)
      } else if (is.matrix(samples)) {
        d <- ncol(samples)
        boundaries <- matrix(NA, nrow = d, ncol = 2)
        for (j in seq_len(d)) {
          min_j <- min(samples[, j])
          max_j <- max(samples[, j])
          padding_j <- 0.1 * (max_j - min_j)
          boundaries[j, ] <- c(min_j - padding_j, max_j + padding_j)
        }
      } else {
        stop("Unsupported samples type: expected vector or matrix.")
      }
    }

    # Check samples and grids dimensionality
    if (is.vector(samples)) {
      dimension <- 1

      sort_index <- order(samples)
      unsort_index <- order(sort_index)  # this gives the positions to map back

      # Sorted samples for computing base measures
      # This is required because further we used trapz for 1 dimensional normalizing
      # constant and trapz function needs ordered input. Instead of doing this
      # inside our functions it's better to do it here.
      # Get the sort index and reverse mapping
      samples <- samples[sort_index]

      if (!is.vector(grids)) stop("If samples is a vector, grids must also be a vector.")
      if (!is.vector(boundaries) || length(boundaries) != 2)
        stop("For 1D input, boundaries must be a vector of length 2 specifying the min and max of the domain.")
      if (length(samples) > length(grids)) {
        warning("grids has fewer points than the samples. Consider increasing grids size for better resolution.")
      }

    } else if (is.matrix(samples)) {
      dimension <- ncol(samples)
      if (!is.matrix(grids) || ncol(grids) != dimension)
        stop("If samples is a matrix, grids must also be a matrix with the same number of columns.")
      if (!is.matrix(boundaries) || nrow(boundaries) != dimension || ncol(boundaries) != 2)
        stop("For d-dimensional input, boundaries must be a d × 2 matrix with each row specifying the min and max of the domain for one dimension.")
      if (nrow(samples) > nrow(grids)) {
        warning("grids has fewer points than the samples. Consider increasing grids size for better resolution.")
      }

    } else {
      stop("samples must be either a numeric vector or a numeric matrix.")
    }

  } else if (data_type == "order") {
    message("Order data detected.")

    # If samples is a list of permutations, convert to numeric matrix
    if (is.list(samples)) {
      samples <- do.call(rbind, samples)  # each element numeric vector
    }

    if (!is.matrix(samples)) {
      samples <- as.matrix(samples)
    }

    mode(samples) <- "numeric"  # ensure numeric matrix

    dimension <- ncol(samples)
  } else if(data_type == "graph"){
    message("Graph data detected: this part of the function has not been implemented yet.")
    # If samples is a list of permutations, convert to numeric matrix
    if (is.list(samples)) {
      samples <- do.call(rbind, samples)  # each element numeric vector
    }

    if (!is.matrix(samples)) {
      samples <- as.matrix(samples)
    }

    mode(samples) <- "numeric"  # ensure numeric matrix

    dimension <- ncol(samples)
  } else{
    stop("Incompatible data type!")
  }


  # Compute the centered kernel matrix at samples
  centered_kernel_mat_samples <- centered_kernel_matrix(
    eval_points_1 = samples,
    eval_points_2 = samples,
    centering_grid = grids,
    hurst_coef = 0.5,
    dimension = dimension,
    data_type = data_type
  )

  #message("centered_kernel_mat_samples is found successfully.")

  # Check if density should be computed on the grids
  density_only_samples <- length(samples) >= 1000
  if (density_only_samples) {
    grids <- samples
    warning("Your samples size is large. To avoid computational complexity,
             the output includes only density estimates at samples points.")
  }

  # Compute kernel matrices for the grids if needed
  if (!density_only_samples) {
    centered_kernel_mat_grids <- centered_kernel_matrix(
      eval_points_1 = samples,
      eval_points_2 = grids,
      centering_grid = grids,
      hurst_coef = 0.5,
      dimension = dimension,
      data_type = data_type
    )

    centered_kernel_self_grids <- diag(centered_kernel_matrix(
      eval_points_1 = grids,
      eval_points_2 = grids,
      centering_grid = grids,
      hurst_coef = 0.5,
      dimension = dimension,
      data_type = data_type
    ))
  }

  # Estimating the base measure
  base_measure_weights <- get_base_measures(samples, data_type = data_type,
                                            dimension = dimension, boundaries)

  #message("base_measure_weights is found successfully.")

  samples_mat <- if (is.matrix(samples)) samples else as.matrix(samples)
  grids_mat   <- if (is.matrix(grids))   grids   else as.matrix(grids)

  # Estimate the weight vector using the Barzilai-Borwein optimization method
  weights_hat <- get_weights(
    lambda = lambda,
    tau = tau,
    centered_kernel_mat_samples = centered_kernel_mat_samples,
    samples = samples,
    base_measure_weights = base_measure_weights,
    dimension = dimension,
    data_type = data_type,
    prior_var_prob = prior_var_prob
  )


  #message("weights_hat is found successfully.")
  #print(weights_hat)
  #print(summary(weights_hat))

  # Compute density estimates based on whether grids evaluation is required
  if (!density_only_samples) {

    base_measure_weights_grids <- get_base_measures(grids, data_type = data_type,
                                              dimension = dimension, boundaries)

    dens <- get_dens(
      centered_kernel_mat_samples,
      centered_kernel_mat_grids,
      centered_kernel_self_grids,
      samples,
      grids,
      base_measure_weights_grids,
      dimension,
      data_type,
      lambda,
      as.vector(weights_hat)
    )

    if(data_type == "euclidean" && dimension == 1){
      # This is required because due to some time saving we sorted samples above
      # and now dens is for sorted index so we need to order it based on input so
      # the user see the correct output.
      dens_samples <- as.vector(dens$samples)[unsort_index]

      dens_unnorm_samples <- as.vector(dens$unnorm_samples)[unsort_index] ##
    }else{
      dens_samples <- dens$samples

      dens_unnorm_samples <- as.vector(dens$unnorm_samples) ##
    }


    # Store results including grids estimates
    result_list <- list(
      weights = as.vector(weights_hat),
      dens_samples = as.vector(dens_samples),
      dens_unnorm_samples = as.vector(dens_unnorm_samples), ##
      dens_grids = as.vector(dens$grids),
      dens_unnorm_grids = as.vector(dens$unnorm_grids),##
      norm_cte = dens$norm_cte,
      data_type = data_type
    )
  } else {
    dens <- get_dens_wo_grid(
      centered_kernel_mat_samples,
      samples = samples,
      base_measure_weights = base_measure_weights,
      dimension = dimension,
      data_type = data_type,
      lambda = lambda,
      weight_vec = as.vector(weights_hat)
    )

    if(data_type == "euclidean" && dimension == 1){
      # This is required because due to some time saving we sorted samples above
      # and now dens is for sorted index so we need to order it based on input so
      # the user see the correct output.
      dens <- dens[unsort_index]
    }

    # Store results (no grids estimates)
    result_list <- list(
      weights = as.vector(weights_hat),
      dens_samples = as.vector(dens)
    )
  }

  #message("density is found successfully.")
  # End timer
  end_time <- Sys.time()

  # Compute total execution time
  total_time <- difftime(end_time, start_time, units = "secs")

  # Store the execution time.
  result_list$time <- total_time

  # Return the result
  return(result_list)
}

