################################################################################
###########################       CLAW  1000       ##############################
################################################################################
library(spatstat)
library(ks)
library(ggplot2)
set.seed(7)

# Define the weights for the mixture distribution
mixture_weights = c(1/2, 1/2)

# Define the parameters for the normal distributions
# First distribution: N(0, 1)

dist_of_peaks <- 10
sd <- 10^-2

means = c(-dist_of_peaks/2,dist_of_peaks/2)
sds = c(sd,sd)

samples <- rnorm_mixture(100, means, sds, mixture_weights)
n <- length(samples)

grids <-  seq(-10.1,10.1,length.out = 4*n)

lambda <- 1
tau <- (lambda^2)/1350

#samples <- sort(samples)
kef_res <- kef(samples,grids = grids,lambda = 1, tau = 1/1350,data_type = "euclidean")

kef_res$time

kef_df <- data.frame(grid = samples, kef_pdf = kef_res$dens_samples)

# Define a matrix of normal densities for each mean and standard deviation
density_matrix <- sapply(seq_along(means), function(i) {
  dnorm(grids, mean = means[i], sd = sds[i])
})

# Define a matrix of normal densities for each mean and standard deviation
density_matrix_samples <- sapply(seq_along(means), function(i) {
  dnorm(samples, mean = means[i], sd = sds[i])
})

# Calculate the true density by taking the weighted sum of the columns
true_density <- density_matrix %*% mixture_weights

# Calculate the true density by taking the weighted sum of the columns
true_density_samples <- density_matrix_samples %*% mixture_weights

true_density_df <- data.frame(grid = grids, true_pdf = true_density)
true_density_df_samples <- data.frame(grid = samples, true_pdf = true_density_samples)

# Perform the adaptive KDE
#kde_adaptive <- akj(sampled_x,sampled_x,kappa = 0.35,alpha = 0.9)
#kde_adaptive <- akj(sampled_x,sampled_x)
#kde_adaptive_df <- data.frame(grid = sampled_x, kde_adaptive_pdf = kde_adaptive$dens)

# Perform the adaptive KDE
kde_adaptive <- densityAdaptiveKernel(samples)

kde_adaptive_df <- data.frame(grid = kde_adaptive$x, kde_adaptive_pdf = kde_adaptive$y)

kde_fixed <- kde(samples,eval.points = samples)
kde_fixed_df <- data.frame(grid = samples, kde_fixed_pdf = kde_fixed$estimate)



ggplot() +
  geom_histogram(aes(x = samples, y = ..density..), fill = 'gray', alpha = 1, color = 'black') +
  geom_line(data = true_density_df_samples, aes(x = grid, y = true_pdf, color = 'True Density'), linewidth = 1) +
  #geom_point(data = true_density_df_sampled, aes(x = grid, y = weights_var, color = 'Weights Var'), size = 1) +
  geom_line(data = kde_adaptive_df, aes(x = grid, y = kde_adaptive_pdf, color = 'KDE Adaptive'), linewidth = 1) +
  geom_line(data = kde_fixed_df, aes(x = grid, y = kde_fixed_pdf, color = 'KDE Fixed'), linewidth = 1) +
  geom_line(data = kef_df, aes(x = grid, y = kef_pdf, color = 'KEF'), linewidth = 1) +
  scale_color_manual(name = "Type of density:", values = c('True Density' = 'red', 'KDE Adaptive' = 'blue','KDE Fixed' = 'limegreen', 'KEF' = 'orange')) +
  #ggtitle(paste('Kernel Density Estimate for lambda_hat =',
  #              format(lambda,digits = 3,scientific = T),'and tau_hat =',format(tau,digits = 3,scientific = T))) +
  xlab('x') +
  ylab('Density') +
  theme_bw() +
  theme(
    legend.position = c(1, 1),             # top-right inside plot
    legend.justification = c("right", "top"),
    legend.background = element_rect(color = alpha("black", 0.6)),
    legend.key.size = unit(1.2, "cm"),           # increase size of legend keys
    legend.text = element_text(size = 14),       # increase text size
    legend.title = element_text(size = 15)       # increase title size
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 3))  # increase line width in legend
  )




