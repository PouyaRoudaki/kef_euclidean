library(dplyr)

lambda_grid <- 10^(seq(-2,2,by=0.2))
tau_grid <- 10^(seq(-8,1,by=0.2))

# Create a data frame with all possible combinations of lambda and tau
grid <- expand.grid(lambda = lambda_grid, tau = tau_grid)

# Calculate log10(tau) and log10(lambda)

grid$log10_lambda <- log10(grid$lambda)
grid$log10_tau <- log10(grid$tau)

# Filter the grid based on the condition
filtered_grid <- grid %>% filter(log10_tau >= log10_lambda - 4.2 )

filtered_grid$mmd_grids <- 0
filtered_grid$ise_grids <- 0
######

# Define the weights for the mixture distribution
mixture_weights = c(1/2, 1/10, 1/10, 1/10, 1/10, 1/10)

# Define the parameters for the normal distributions
# First distribution: N(0, 1)
means = c(0, -1, -0.5, 0, 0.5, 1)
sds = c(1, 0.1, 0.1, 0.1, 0.1, 0.1)

samples <- rnorm_mixture(100, means, sds, mixture_weights)
n <- length(samples)

grids <-  seq(-3.1,3.1,length.out = 4*n)

# Mixture Normal Example
# Define a matrix of normal densities for each mean and standard deviation
density_matrix_grids <- sapply(seq_along(means), function(i) {
  dnorm(grids, mean = means[i], sd = sds[i])
})

# Calculate the true density by taking the weighted sum of the columns
true_density_grids <- as.numeric(density_matrix_grids %*% mixture_weights)

true_density_df_grids <- data.frame(grid = grids, true_pdf = true_density_grids)

centered_kernel_self_grids_mat <- centered_kernel_matrix(dimension = 1, grids, grids, grids, 0.5)

#####


for (i in 462:nrow(filtered_grid)) {

  lambda <- filtered_grid$lambda[i]
  tau <- filtered_grid$tau[i]

  kef_res <- kef(samples,grids = grids,lambda,tau)

  estimated_density_grids <- kef_res$dens_grids

  ise_grids <- l2_ise(true_density_grids, estimated_density_grids)
  mmd_grids <- mmd_grids(centered_kernel_self_grids_mat, true_density_grids, estimated_density_grids)

  filtered_grid$ise_grids[i] <- ise_grids
  filtered_grid$mmd_grids[i] <- mmd_grids

}

#filtered_grid[filtered_grid$mmd_grids == 0 | filtered_grid$ise_grids == 0,c(5,6)] <- 10000


best_mmd_each_lambda <- filtered_grid %>%
  group_by(lambda) %>%
  slice_min(mmd_grids, with_ties = FALSE)


best_mmd_each_lambda <- filtered_grid %>%
  group_by(lambda) %>%
  slice_min(ise_grids, with_ties = FALSE)
