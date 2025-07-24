library(MASS)
library(ks)
library(ggplot2)

set.seed(7)

# === Generate samples from a standard 2D normal ===
n <- 100
mu <- c(0, 0)
Sigma <- matrix(c(1, 0.8, 0.8, 1), 2, 2)  # Correlated normal
samples <- mvrnorm(n = n, mu = mu, Sigma = Sigma)

####
set.seed(7)  # for reproducibility

# Number of total samples
n_total <- 100

# Mixture proportions (e.g., 60% from first, 40% from second)
p1 <- 0.6
p2 <- 0.4

# Sample sizes for each component
n1 <- round(n_total * p1)
n2 <- n_total - n1

# --- First component ---
mu1 <- c(1, -2)
Sigma1 <- matrix(c(1, 0.8, 0.8, 1), 2, 2)
samples1 <- mvrnorm(n = n1, mu = mu1, Sigma = Sigma1)

# --- Second component ---
mu2 <- c(-1, 2)
Sigma2 <- matrix(c(1, -0.7, -0.7, 1), 2, 2)
samples2 <- mvrnorm(n = n2, mu = mu2, Sigma = Sigma2)

# Combine into mixture
samples <- rbind(samples1, samples2)

# Optional: Plot
plot(samples, col = rep(c("blue", "red"), times = c(n1, n2)),
     pch = 19, xlab = "X1", ylab = "X2", main = "Mixture of Two 2D Normals")
legend("topright", legend = c("Component 1", "Component 2"), col = c("blue", "red"), pch = 19)

# === Fixed KDE ===
H <- Hpi(samples)  # bandwidth selection
kde_fixed <- kde(x = samples, H = H, compute.cont = TRUE)

# === Create grid and evaluate KDE ===
x_grid <- seq(-4.5, 4.5, length.out = 800)
y_grid <- seq(-4.5, 4.5, length.out = 800)
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

######
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

#install.packages("plot3D")
library(plot3D)


# Example data (if you don’t already have them defined)
# samples <- matrix(runif(200), ncol = 2)
# sample1_kernel <- runif(nrow(samples))

# Extract coordinates and values
x <- samples[,1]
y <- samples[,2]
z1 <- sample1_kernel
z2 <- sample2_kernel

dx <- dy <- 0.1

# Set plot limits
xlim <- range(x) + c(-dx, dx)
ylim <- range(y) + c(-dy, dy)
zlim <- c(0, max(z1) * 1.1)

# Initialize empty 3D plot
scatter3D(x=0, y=0, z=0, xlim = xlim, ylim = ylim, zlim = zlim,
          pch = 1, cex = 0.1, col = "black",
          xlab = "X", ylab = "Y", zlab = "Z",
          main = "3D Barplot of Sample Points",
          theta = 45,  # rotate horizontally
          phi = 30 )

# Add bars
for (i in 1:length(x)) {
  box3D(x0 = x[i] - dx/2, y0 = y[i] - dy/2, z0 = 0,
         x1 = x[i] + dx/2, y1 = y[i] + dy/2, z1 = z1[i],
         col = "skyblue", border = "black", add = TRUE)
}

points3D(x = x[1], y = y[1], z = 0,
         col = "red",       # red fill color
         pch = 19,          # solid round point
         cex = 3,         # smaller size
         add = TRUE)

# Initialize empty 3D plot
scatter3D(x=0, y=0, z=0, xlim = xlim, ylim = ylim, zlim = zlim,
          pch = 1, cex = 0.1, col = "black",
          xlab = "X", ylab = "Y", zlab = "Z",
          main = "3D Barplot of Sample Points",
          theta = 45,  # rotate horizontally
          phi = 30 )

# Add bars
for (i in 1:length(x)) {
  box3D(x0 = x[i] - dx/2, y0 = y[i] - dy/2, z0 = 0,
        x1 = x[i] + dx/2, y1 = y[i] + dy/2, z1 = z2[i],
        col = "skyblue", border = "black", add = TRUE)
}

points3D(x = x[2], y = y[2], z = 0,
         col = "red",       # red fill color
         pch = 19,          # solid round point
         cex = 3,         # smaller size
         add = TRUE)

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

d <- ncol(samples)
boundaries <- matrix(NA, nrow = d, ncol = 2)
for (j in seq_len(d)) {
  min_j <- min(samples[, j])
  max_j <- max(samples[, j])
  padding_j <- 0.1 * (max_j - min_j)
  boundaries[j, ] <- c(min_j - padding_j, max_j + padding_j)
}
print(boundaries)

# Estimating the base measure
base_measure_weights <- get_base_measures(samples, boundaries, dimension = dimension)

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


#get_dens_or_prob(centered_kernel_mat_samples,
#                 centered_kernel_mat_grids,
#                 centered_kernel_self_grids,
#                 samples,
#                 grids,
#                 lambda_hat,
#                 type_of_p_is_prob = TRUE,
#                 type_of_q_is_prob = TRUE,
#                 weight_hat_vec)



get_dens_or_prob <- function(centered_kernel_mat_samples,
                               centered_kernel_mat_grids,
                               centered_kernel_self_grids,
                               samples,
                               grids,
                               lambda_hat,
                               type_of_p_is_prob = TRUE,
                               type_of_q_is_prob = TRUE,
                               weight_hat_vec){

  # Compute the probabilities at the sampled points
  dens_samples <- density_at_samples(centered_kernel_mat_samples, lambda_hat, weight_hat_vec)

  print(dim(dens_samples))
  print(length(dens_samples))
  print("***")
  # Compute the densities at the grid points
  dens_grids <- density_at_grids(centered_kernel_mat_grids, centered_kernel_self_grids, lambda_hat, weight_hat_vec)

  print(dim(dens_grids))
  print(length(dens_grids))
  print("***")

  # Normalize the density by the integral over the grid
  normalizing_cte <- trapz2d( grids[,1], grids[,2], dens_grids)  # trapz is from pracma package

  # Prepare the output as a list of normalized probabilities
  dens_list <- list()
  dens_list$samples <- dens_samples / normalizing_cte
  dens_list$grids <- dens_grids / normalizing_cte

  prob_list <- list()
  prob_list$samples <- dens_list$samples / sum(dens_list$samples)
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
                              grids,
                              lambda_hat,
                              type_of_p_is_prob = TRUE,
                              type_of_q_is_prob = TRUE,
                              weight_hat_vec)

    prob_samples <- probs$samples
    prob_grids <- probs$grids

    s <- lambda_hat * (colSums(centered_kernel_mat_samples) -
                         n * prob_grids %*% t(centered_kernel_mat_grids)) -
      tau_hat * weight_hat_vec / prob_samples

    # Compute the inverse of the Hessian matrix

    Hessian <- lambda_hat^2 * n * (centered_kernel_mat_grids %*%
                                   diag(prob_grids) %*% t(centered_kernel_mat_grids) -
                                   (centered_kernel_mat_grids %*% prob_grids) %*%
                                  t(centered_kernel_mat_grids %*% prob_grids)  ) + diag(tau_hat / prob_samples)

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


######################################################
# Estimate the weight vector using the Barzilai-Borwein optimization method
#weights_hat <- get_weights(
#  lambda = 1,
#  tau = 10,
#  centered_kernel_mat_samples = centered_kernel_mat_samples,
#  samples = samples,
#  base_measure_weights = base_measure_weights,
#  dimension = dimension
#)

#as.numeric(weight_hat_vec) -weights_hat
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

####################

x_grid <- seq(-5, 5, length.out = 100)
y_grid <- seq(-5, 5, length.out = 100)
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

density2d_estimate <- get_dens_or_prob(centered_kernel_mat_samples,centered_kernel_mat_grids,
                 centered_kernel_self_grids,samples,grids, lambda_hat = 1,
                 type_of_p_is_prob = FALSE,
                 type_of_q_is_prob = FALSE,weight_hat_vec = as.numeric(weight_hat_vec))


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


library(ggplot2)
library(scales)  # for scientific()

# Base plot with first KDE contour and points
ggplot() +
  # First contour (black lines from kde_df)
  geom_contour(data = kde_df, aes(x = x, y = y, z = density),
               color = "black", linewidth = 0.5) +

  # Sample points
  geom_point(data = samples_df, aes(x = x, y = y),
             color = "black", alpha = 0.3, size = 0.8) +

  # Second contour (colored lines from df)
  geom_contour(data = df, aes(x = x, y = y, z = weight, color = ..level..),
               bins = 15) +

  # Color scale for second contour
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
    title = "KEF estimation lambda = 1, tau = 0.01",
    #title = "KDE estimation",
    subtitle = paste("n = 100,", nrow(df), "grid points"),
    x = "X",
    y = "Y"
  ) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

