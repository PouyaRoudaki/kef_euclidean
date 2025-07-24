library(MASS)
library(ks)
library(ggplot2)

set.seed(7)

# === Generate samples from a standard 2D normal ===
n <- 100
mu <- c(0, 0)
Sigma <- matrix(c(1, 0.8, 0.8, 1), 2, 2)  # Correlated normal
samples <- mvrnorm(n = n, mu = mu, Sigma = Sigma)

# === Fixed KDE ===
H <- Hpi(samples)  # bandwidth selection
kde_fixed <- kde(x = samples, H = H, compute.cont = TRUE)

# === Create grid and evaluate KDE ===
x_grid <- seq(-4, 4, length.out = 30)
y_grid <- seq(-4, 4, length.out = 30)
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
x_grid <- seq(-2.5, 2.5, length.out = 40)
y_grid <- seq(-2.5, 2.5, length.out = 40)
grid_mat <- as.matrix(expand.grid(x_grid, y_grid))

# KEF estimation
kef_res <- kef(samples = samples, grids = grid_mat, lambda = 1, tau = 0.01)


#library(scatterplot3d)

#scatterplot3d(
#  x = samples[,1],
#  y = samples[,2],
#  z = kef_res$weights,
#  pch = 16,
#  color = "blue",
#  main = "3D Scatterplot of Weights",
#  xlab = "X1",
#  ylab = "X2",
#  zlab = "Weight"
#)

#library(rgl)

#plot3d(
#  x = samples[,1],
#  y = samples[,2],
#  z = kef_res$weights,
#  col = "red",
#  size = 1,
#  type = "s",  # "s" for spheres, "p" for points
#  xlab = "X1",
#  ylab = "X2",
#  zlab = "Weight"
#)

#plot3d(
#  x = samples[,1],
#  y = samples[,2],
#  z = kef_res$dens_samples,
#  col = "red",
#  size = 1,
#  type = "s",  # "s" for spheres, "p" for points
#  xlab = "X1",
#  ylab = "X2",
#  zlab = "Dens"
#)


kef_df <- data.frame(
  x = grid_mat[,1],
  y = grid_mat[,2],
  density = kef_res$dens_grids
)


# Scatter + contour plot
ggplot(kef_df, aes(x = x, y = y)) +
  geom_contour_filled(aes(z = density), color = "black", linewidth = 1) +
  geom_point(data = samples_df, aes(x = x, y = y), color = "black", alpha = 0.6, size = 1) +
  scale_fill_viridis_d() +
  labs(title = "2D KEF Density Contours", x = "X1", y = "X2") +
  theme_bw()



