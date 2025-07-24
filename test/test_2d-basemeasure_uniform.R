library(MASS)
library(ks)
library(ggplot2)

set.seed(7)

# === Generate samples from a 2D uniform ===
n <- 100

samples <- cbind(runif(n), runif(n))

plot(samples, col = rep(c("blue"), times = c(n)),
     pch = 19, xlab = "X1", ylab = "X2", main = "2D unifrom samples")


H <- Hpi(samples)  # bandwidth selection
kde_fixed <- kde(x = samples, H = H, compute.cont = TRUE)

# === Create grid and evaluate KDE ===
x_grid <- seq(0, 1, length.out = 800)
y_grid <- seq(0, 1, length.out = 800)
grid_df <- expand.grid(x = x_grid, y = y_grid)

# Evaluate KDE at grid points
z_kde <- predict(kde_fixed, x = grid_df)

# Combine into data frame for ggplot
kde_df <- cbind(grid_df, density = as.vector(z_kde))

# === Plot KDE as contour + sample points ===
# Convert samples to data frame with named columns
samples_df <- data.frame(x = samples[, 1], y = samples[, 2])

# Plot
ggplot(kde_df, aes(x = x, y = y)) +
  geom_contour_filled(aes(z = density), alpha = 0.8) +
  geom_point(data = samples_df, aes(x = x, y = y), color = "black", alpha = 0.3, size = 0.8) +
  scale_fill_viridis_d() +
  labs(title = "2D KDE from standard normal", fill = "Density") +
  theme_minimal()

# Create contour plot with lines only
ggplot(kde_df, aes(x = x, y = y)) +
  geom_contour(aes(z = density), color = "black", linewidth = 1) +
  geom_point(data = samples_df, aes(x = x, y = y), color = "black", alpha = 0.3, size = 0.8) +
  labs(title = "2D KDE from standard normal") +
  theme_bw()

#####


# === Create grid and evaluate KDE ===
x_grid <- seq(-5, 5, length.out = 50)
y_grid <- seq(-5, 5, length.out = 50)
grid_mat <- as.matrix(expand.grid(x_grid, y_grid))

dimension <- 2

grids <- grid_mat

# Compute the centered kernel matrix at samples
centered_kernel_mat_samples <- centered_kernel_matrix(
  dimension = dimension,
  eval_points_1 = samples,
  eval_points_2 = samples,
  centering_grid = grids,
  hurst_coef = 0.5
)

sample1_kernel <- centered_kernel_mat_samples[1,]
sample2_kernel <- centered_kernel_mat_samples[2,]


d <- ncol(samples)
boundaries <- matrix(c(0,0,1,1),nrow = 2)

print(boundaries)

# Estimating the base measure
base_measure_weights <- get_base_measures(samples, boundaries, dimension = dimension)

print(base_measure_weights)

### Check base measure
x1 <- samples[, 1]
x2 <- samples[, 2]

# Compute Voronoi diagram
vor <- deldir(x1, x2, rw = c(0, 1, 0, 1))
tiles <- tile.list(vor)
summary_df <- vor$summary

# Start empty plot
plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1), asp = 1,
     xlab = "x1", ylab = "x2", main = "Voronoi Diagram with Area Labels")

# Draw all tiles
for (i in seq_along(tiles)) {
  tile <- tiles[[i]]

  # Color the first sample’s tile green
  if (i == 1) {
    polygon(tile$x, tile$y, col = "lightgreen", border = "black")
  } else {
    polygon(tile$x, tile$y, border = "black")
  }

  # Centroid for area label
  cx <- mean(tile$x)
  cy <- mean(tile$y)
  area_label <- round(summary_df$dir.area[i], 4)
  text(cx, cy, labels = area_label, cex = 0.7, col = "red")
}

# Add the sample points
points(x1, x2, pch = 16, col = "blue")

# Optionally highlight sample[1,] point with larger green point
points(x1[1], x2[1], pch = 16, col = "darkgreen", cex = 1.5)



##################################################################


# Get density
density_at_samples <- function(centered_kernel_mat_samples, lambda, weight_hat_vec) {

  # Extract the diagonal elements of the kernel matrix, which represent the
  # self-kernel values (i.e., the kernel values of each point with itself).
  diag_vals <- diag(centered_kernel_mat_samples)

  # Compute the density values at each sampled data point
  p_vec <- sapply(1:nrow(centered_kernel_mat_samples), function(i) {
    # For each sampled point, calculate the exponential of the linear combination of weights
    # and centered kernel values, adjusted by lambda_hat, and subtract 0.5 times the
    # self-kernel value for that point.
    exp(lambda_hat * (weight_hat_vec %*% centered_kernel_mat_samples[,i] -
                        0.5 * diag_vals[i]))
  })

  # Return the vector of density values
  return(p_vec)
}

density_at_grids <- function(centered_kernel_mat_grids,
                             centered_kernel_self_grids,
                             lambda_hat, weight_hat_vec) {

  # Compute the density values at each grid point
  p_vec <- sapply(1:ncol(centered_kernel_mat_grids), function(i) {
    # For each grid point, calculate the exponential of the linear combination of weights
    # and centered kernel values, adjusted by lambda_hat, and subtract 0.5 times the
    # self-grid kernel value.
    exp(lambda_hat * (weight_hat_vec %*% centered_kernel_mat_grids[,i] -
                        0.5 * centered_kernel_self_grids[i]))
  })

  # Return the vector of density values
  return(p_vec)
}

trapz2d <- function(x, y, z, byrow = FALSE) {
  # x: vector of length n (x-axis)
  # y: vector of length m (y-axis)
  # z: vector of length m * n (function values)
  # byrow: whether z is ordered row-wise (default FALSE = column-wise)
  x <- unique(x)
  y <- unique(y)

  n <- length(x)
  m <- length(y)

  if (length(z) != m * n) {
    stop("Length of z must be equal to length(x) * length(y).")
  }

  # Reshape z to matrix: m rows (y), n columns (x)
  z_mat <- matrix(z, nrow = m, ncol = n, byrow = byrow)

  # Now apply the 2D trapezoidal rule
  dx <- diff(x)
  dy <- diff(y)
  integral <- 0

  for (i in 1:(m - 1)) {
    for (j in 1:(n - 1)) {
      # Get average over the cell's 4 corners
      f_avg <- (z_mat[i, j] + z_mat[i+1, j] + z_mat[i, j+1] + z_mat[i+1, j+1]) / 4
      area <- dx[j] * dy[i]
      integral <- integral + f_avg * area
    }
  }

  return(integral)
}

get_dens_or_prob <- function(centered_kernel_mat_samples,
                             centered_kernel_mat_grids,
                             centered_kernel_self_grids,
                             samples,
                             base_measure_weights,
                             lambda_hat,
                             type_of_p_is_prob = TRUE,
                             weight_hat_vec){

  # Compute the probabilities at the sampled points
  dens_samples <- density_at_samples(centered_kernel_mat_samples, lambda_hat, weight_hat_vec)

  #print(dim(dens_samples))
  print(length(dens_samples))
  print("***")
  # Compute the densities at the grid points
  #dens_grids <- density_at_grids(centered_kernel_mat_grids, centered_kernel_self_grids, lambda_hat, weight_hat_vec)

  #print(dim(dens_grids))
  #print(length(dens_grids))
  #print("***")

  # Normalize the density by the integral over the grid
  #normalizing_cte <- trapz2d( grids[,1], grids[,2], dens_grids)  # trapz is from pracma package
  normalizing_cte <- dens_samples %*% base_measure_weights

  print(dim(normalizing_cte))
  print(length(normalizing_cte))

  normalizing_cte <- as.numeric(normalizing_cte)
  # Prepare the output as a list of normalized probabilities
  dens_list <- list()
  dens_list$samples <- dens_samples / normalizing_cte
  #dens_list$grids <- dens_grids / normalizing_cte

  prob_list <- list()
  prob_list$samples <- (base_measure_weights * dens_list$samples) / sum((base_measure_weights * dens_list$samples))
  #prob_list$grids <- dens_list$grids / sum(dens_list$grids)

  #print(paste("inside get_dens_or_probs",length(prob_list$grid_x)))
  #if(method_of_p_calculation == "neighborhood_grid"){

  #  approx_dens_or_prob <- get_grid_approx_dens_or_probs_vectorized(sampled_x, x_grid, dens_list)

  #  dens_list$sampled_x <- approx_dens_or_prob$dens

  #  prob_list$sampled_x <- approx_dens_or_prob$prob

  #}

  result_list <- list()

  if(type_of_p_is_prob){
    result_list$samples <- prob_list$samples
  }else{
    result_list$samples <- dens_list$samples
  }

  #if(type_of_q_is_prob){
  #  result_list$grids <- prob_list$grids
  #}else{
  #  result_list$grids <- dens_list$grids
  #}

  #print(paste("inside get_dens_or_probs",length(result_list$grid_x)))

  return(result_list)
}

# Get weights
lambda_hat <- 1
tau_hat <- 0.01
print_trace = TRUE

max_iteration <- 300  # Maximum number of iterations for the Newton-Raphson method
NRstepsize <- 0.1  # Step size for the Newton-Raphson update
n <- nrow(centered_kernel_mat_samples)  # Number of sampled points
weight_hat_vec <- rep(0, n)  # Initialize the weight vector with zeros
#s <- rep(1000, n)
#weight_hat_change <- rep(1000, n)
#counter <- 1
#while ((norm(s, p = 2) > 10^(-10)) & (norm(weight_hat_change, p = 2) > 10^(-10))) {


for (i in 1:max_iteration) {
  # Calculate probabilities for sampled and grid points
  probs <- get_dens_or_prob(centered_kernel_mat_samples,
                            centered_kernel_mat_grids,
                            centered_kernel_self_grids,
                            samples,
                            base_measure_weights,
                            lambda_hat,
                            type_of_p_is_prob = TRUE,
                            weight_hat_vec)

  prob_samples <- probs$samples

  s <- lambda_hat * (colSums(centered_kernel_mat_samples) -
                       n * prob_samples %*% t(centered_kernel_mat_samples)) -
    tau_hat * weight_hat_vec / prob_samples

  # Compute the inverse of the Hessian matrix

  Hessian <- lambda_hat^2 * n * (centered_kernel_mat_samples %*%
                                   diag(prob_samples) %*% t(centered_kernel_mat_samples) -
                                   (centered_kernel_mat_samples %*% prob_samples) %*%
                                   t(centered_kernel_mat_samples %*% prob_samples)  ) + diag(tau_hat / prob_samples)

  Hessian_inv <- solve(Hessian)


  #print(s)
  #print(Hessian_inv)
  # Update the weight vector using the Newton-Raphson method
  #weight_hat_change <- NRstepsize * s %*% Hessian_inv
  weight_hat_vec <- weight_hat_vec + NRstepsize * s %*% Hessian_inv

  # Print progress every 10% of the iterations or at the first iteration
  #if ((i %% round(max_iteration / 10) == 0 || i == 1) & print_trace == TRUE) {
  print(paste("Iteration", i, ": ||s||_2 =", pracma::Norm(s)))
  #print(summary(as.vector(Hessian)))
  #}
  #counter = counter + 1
}


s.numeric(weight_hat_vec) -weights_hat
library(ggplot2)
library(scales)
library(viridis)

# Example data (replace with your own)
# samples <- matrix(c(...), ncol = 2)
# weights_hat <- c(...)

# Create data frame
#df <- data.frame(
#  x = samples[, 1],
#  y = samples[, 2],
#  weight = weights_hat
#)

df <- data.frame(
  x = samples[, 1],
  y = samples[, 2],
  weight = as.numeric(weight_hat_vec)
)

# Basic weighted scatter plot with red-blue spectrum
p1 <- ggplot(df, aes(x = x, y = y, size = abs(weight), color = weight)) +
  geom_point(alpha = 0.8) +
  scale_color_gradient2(
    low = "blue",
    mid = "purple",
    high = "red",
    midpoint = -0.1,
    name = "Weight",
    labels = scientific
  ) +
  scale_size_continuous(
    range = c(2, 15),
    name = "Weight",
    labels = scientific
  ) +
  theme_minimal() +
  labs(
    title = "Weighted Scatter Plot",
    subtitle = paste("n =", nrow(df), "points"),
    x = "X",
    y = "Y"
  ) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )+
  theme_bw()

print(p1)

##################################################

get_dens_or_prob_just_for_output <- function(centered_kernel_mat_samples,
                             centered_kernel_mat_grids,
                             centered_kernel_self_grids,
                             samples,
                             base_measure_weights,
                             grids,
                             lambda_hat,
                             type_of_p_is_prob = TRUE,
                             type_of_q_is_prob = TRUE,
                             weight_hat_vec){

  # Compute the probabilities at the sampled points
  dens_samples <- density_at_samples(centered_kernel_mat_samples, lambda_hat, weight_hat_vec)

  #print(dim(dens_samples))
  print(length(dens_samples))
  print("***")
  # Compute the densities at the grid points
  dens_grids <- density_at_grids(centered_kernel_mat_grids, centered_kernel_self_grids, lambda_hat, weight_hat_vec)

  #print(dim(dens_grids))
  print(length(dens_grids))
  print("***")

  # Normalize the density by the integral over the grid
  normalizing_cte <- trapz2d( grids[,1], grids[,2], dens_grids)  # trapz is from pracma package
  #normalizing_cte <- dens_samples %*% base_measure_weights

  print(dim(normalizing_cte))
  print(length(normalizing_cte))

  normalizing_cte <- as.numeric(normalizing_cte)
  # Prepare the output as a list of normalized probabilities
  dens_list <- list()
  dens_list$samples <- dens_samples / normalizing_cte
  dens_list$grids <- dens_grids / normalizing_cte

  prob_list <- list()
  prob_list$samples <- (base_measure_weights * dens_list$samples) / sum((base_measure_weights * dens_list$samples))
  prob_list$grids <- dens_list$grids / sum(dens_list$grids)

  #print(paste("inside get_dens_or_probs",length(prob_list$grid_x)))
  #if(method_of_p_calculation == "neighborhood_grid"){

  #  approx_dens_or_prob <- get_grid_approx_dens_or_probs_vectorized(sampled_x, x_grid, dens_list)

  #  dens_list$sampled_x <- approx_dens_or_prob$dens

  #  prob_list$sampled_x <- approx_dens_or_prob$prob

  #}

  result_list <- list()

  if(type_of_p_is_prob){
    result_list$samples <- prob_list$samples
  }else{
    result_list$samples <- dens_list$samples
  }

  if(type_of_q_is_prob){
    result_list$grids <- prob_list$grids
  }else{
    result_list$grids <- dens_list$grids
  }

  #print(paste("inside get_dens_or_probs",length(result_list$grid_x)))

  return(result_list)
}

x_grid <- seq(0, 1, length.out = 100)
y_grid <- seq(0, 1, length.out = 100)
grid_mat <- as.matrix(expand.grid(x_grid, y_grid))

dimension <- 2

grids <- grid_mat

centered_kernel_mat_grids <- centered_kernel_matrix(
  dimension = dimension,
  eval_points_1 = samples,
  eval_points_2 = grids,
  centering_grid = grids,
  hurst_coef = 0.5
)

centered_kernel_self_grids <- diag(centered_kernel_matrix(
  dimension = dimension,
  eval_points_1 = grids,
  eval_points_2 = grids,
  centering_grid = grids,
  hurst_coef = 0.5
))

density2d_estimate <- get_dens_or_prob_just_for_output(centered_kernel_mat_samples,
                                                       centered_kernel_mat_grids,
                                       centered_kernel_self_grids,
                                       samples,
                                       base_measure_weights,
                                       grids,
                                       lambda_hat = 1,
                                       type_of_p_is_prob = FALSE,
                                       type_of_q_is_prob = FALSE,
                                       weight_hat_vec = as.numeric(weight_hat_vec))


df <- data.frame(
  x = grids[, 1],
  y = grids[, 2],
  weight = density2d_estimate$grids
)

# Basic weighted scatter plot with red-blue spectrum
p1 <- ggplot(df, aes(x = x, y = y, size = abs(weight), color = weight)) +
  geom_point(alpha = 0.1) +
  scale_color_gradient2(
    low = "blue",
    mid = "purple",
    high = "red",
    midpoint = -0.1,
    name = "Weight",
    labels = scientific
  ) +
  scale_size_continuous(
    range = c(2, 15),
    name = "Weight",
    labels = scientific
  ) +
  theme_minimal() +
  labs(
    title = "Weighted Scatter Plot",
    subtitle = paste("n =", nrow(df), "points"),
    x = "X",
    y = "Y"
  ) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )+
  theme_bw()

print(p1)

# Use contour plot
p_contour <- ggplot(df, aes(x = x, y = y, z = weight)) +
  geom_contour(aes(color = ..level..), bins = 15) +
  scale_color_gradient2(
    low = "blue",
    mid = "purple",
    high = "red",
    midpoint = -0.1,
    name = "Density",
    labels = scientific
  ) +
  coord_fixed() +
  theme_bw() +
  labs(
    title = "Contour Plot of Density Estimate",
    subtitle = paste("n =100,", nrow(df), "grid points, weight grid points = 2500"),
    x = "X",
    y = "Y"
  ) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

print(p_contour)

