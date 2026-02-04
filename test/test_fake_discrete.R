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

dist_of_peaks <- 1000
sd <- 1

means = c(-dist_of_peaks/2,dist_of_peaks/2)
sds = c(sd,sd)

samples <- rnorm_mixture(100, means, sds, mixture_weights)
n <- length(samples)

grids <-  seq(-dist_of_peaks,dist_of_peaks,length.out = 4*n)

lambda <- 1
tau <- (lambda^2)/1350

#samples <- sort(samples)
kef_res <- kef(samples,grids = grids,lambda = lambda, tau = tau,data_type = "euclidean")

neg_idx <- which(samples < 0)
neg_idx[ which.max(samples[neg_idx]) ]
samples[neg_idx[ which.max(samples[neg_idx]) ]]

pos_idx <- which(samples > 0)
pos_idx[ which.min(samples[pos_idx]) ]
samples[pos_idx[ which.min(samples[pos_idx]) ]]




kef_res$time
################################################################################

library(ggplot2)
library(scales)
library(grid)

## 1) plot(samples, kef_res$dens_samples)

df_dens_samples <- data.frame(
  x   = samples,
  y   = kef_res$dens_samples
)

ggplot(df_dens_samples_ord, aes(x = x, y = y)) +
  geom_point(color = "orange", size = 1) +
  geom_line(color = "orange", linewidth = 1) +
  xlab("x") +
  ylab("Estimated density at samples") +
  theme_bw() +
  theme(
    legend.position      = "none",
    legend.background    = element_rect(color = alpha("black", 0.6)),
    legend.key.size      = unit(1.2, "cm"),
    legend.text          = element_text(size = 14),
    legend.title         = element_text(size = 15)
  )

###############################################################################
## 2) plot(samples, kef_res$weights)

df_weights <- data.frame(
  x       = samples,
  weight  = kef_res$weights
)

neg_idx <- which(samples < 0)
i_neg   <- neg_idx[ which.max(samples[neg_idx]) ]

pos_idx <- which(samples > 0)
i_pos   <- pos_idx[ which.min(samples[pos_idx]) ]

df_weights$color <- "blue"
df_weights$color[c(i_neg, i_pos)] <- "red"

df_weights$size <- 2
df_weights$size[c(i_neg, i_pos)] <- 2

df_weights$alpha <- 0.1
df_weights$alpha[c(i_neg, i_pos)] <- 1




idx_two_smallest <- order(df_weights$weight)[1:2]
samples[idx_two_smallest]

sorted <- sort(samples)
which(sorted %in% samples[idx_two_smallest])

min(samples)
max(samples)

ggplot(df_weights, aes(x = x, y = weight)) +
  geom_point(aes(color = color, size = size, alpha = alpha)) +
  scale_color_identity() +     # tells ggplot to use the color values directly
  xlab("x") +
  ylab("KEF weights") +
  theme_bw() +
  theme(
    legend.position      = "none",
    legend.background    = element_rect(color = alpha("black", 0.6)),
    legend.key.size      = unit(1.2, "cm"),
    legend.text          = element_text(size = 14),
    legend.title         = element_text(size = 15)
  )

###############################################################################
df_dens_grids <- data.frame(
  grid = grids,
  dens = kef_res$dens_grids
)



ggplot(df_dens_grids, aes(x = grid, y = dens)) +
  geom_line(color = "orange", linewidth = 1) +
  xlab("x") +
  ylab("Estimated density on grid") +
  theme_bw() +
  theme(
    legend.position      = "none",
    legend.background    = element_rect(color = alpha("black", 0.6)),
    legend.key.size      = unit(1.2, "cm"),
    legend.text          = element_text(size = 14),
    legend.title         = element_text(size = 15)
  )
################################################################################


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

dist_of_peaks <- 100
sd <- 0.01

means = c(-dist_of_peaks/2,dist_of_peaks/2)
sds = c(sd,sd)

samples <- rnorm_mixture(100, means, sds, mixture_weights)
n <- length(samples)

grids <-  seq(-dist_of_peaks,dist_of_peaks,length.out = 4*n)

lambda <- 1
tau <- (lambda^2)/1350

#samples <- sort(samples)
kef_res <- kef(samples,grids = grids,lambda = lambda, tau = tau,data_type = "euclidean")

neg_idx <- which(samples < 0)
neg_idx[ which.max(samples[neg_idx]) ]
samples[neg_idx[ which.max(samples[neg_idx]) ]]

pos_idx <- which(samples > 0)
pos_idx[ which.min(samples[pos_idx]) ]
samples[pos_idx[ which.min(samples[pos_idx]) ]]




kef_res$time
################################################################################

library(ggplot2)
library(scales)
library(grid)

## 1) plot(samples, kef_res$dens_samples)

df_dens_samples <- data.frame(
  x   = samples,
  y   = kef_res$dens_samples
)

ggplot(df_dens_samples_ord, aes(x = x, y = y)) +
  geom_point(color = "orange", size = 1) +
  geom_line(color = "orange", linewidth = 1) +
  xlab("x") +
  ylab("Estimated density at samples") +
  theme_bw() +
  theme(
    legend.position      = "none",
    legend.background    = element_rect(color = alpha("black", 0.6)),
    legend.key.size      = unit(1.2, "cm"),
    legend.text          = element_text(size = 14),
    legend.title         = element_text(size = 15)
  )

###############################################################################
## 2) plot(samples, kef_res$weights)

df_weights <- data.frame(
  x       = samples,
  weight  = kef_res$weights,
  base_measure = get_base_measures(samples,boundaries = c(min(grids), max(grids)))
)

neg_idx <- which(samples < 0)
i_neg   <- neg_idx[ which.max(samples[neg_idx]) ]

pos_idx <- which(samples > 0)
i_pos   <- pos_idx[ which.min(samples[pos_idx]) ]

df_weights$color <- "blue"
df_weights$color[c(i_neg, i_pos)] <- "red"

df_weights$size <- 2
df_weights$size[c(i_neg, i_pos)] <- 2

df_weights$alpha <- 0.1
df_weights$alpha[c(i_neg, i_pos)] <- 1




idx_two_smallest <- order(df_weights$weight)[1:2]
samples[idx_two_smallest]

sorted <- sort(samples)
which(sorted %in% samples[idx_two_smallest])

min(samples)
max(samples)

ggplot(df_weights, aes(x = x, y = weight)) +
  geom_point(aes(color = color, size = size, alpha = alpha)) +
  scale_color_identity() +     # tells ggplot to use the color values directly
  xlab("x") +
  ylab("KEF weights") +
  theme_bw() +
  theme(
    legend.position      = "none",
    legend.background    = element_rect(color = alpha("black", 0.6)),
    legend.key.size      = unit(1.2, "cm"),
    legend.text          = element_text(size = 14),
    legend.title         = element_text(size = 15)
  )


ggplot(df_weights, aes(x = x, y = base_measure)) +
  geom_point(aes(color = color, size = size, alpha = alpha)) +
  scale_color_identity() +     # tells ggplot to use the color values directly
  xlab("x") +
  ylab("KEF weights") +
  theme_bw() +
  theme(
    legend.position      = "none",
    legend.background    = element_rect(color = alpha("black", 0.6)),
    legend.key.size      = unit(1.2, "cm"),
    legend.text          = element_text(size = 14),
    legend.title         = element_text(size = 15)
  )

###############################################################################
df_dens_grids <- data.frame(
  grid = grids,
  dens = kef_res$dens_grids
)



ggplot(df_dens_grids, aes(x = grid, y = dens)) +
  geom_line(color = "orange", linewidth = 1) +
  xlab("x") +
  ylab("Estimated density on grid") +
  theme_bw() +
  theme(
    legend.position      = "none",
    legend.background    = element_rect(color = alpha("black", 0.6)),
    legend.key.size      = unit(1.2, "cm"),
    legend.text          = element_text(size = 14),
    legend.title         = element_text(size = 15)
  )


################################################################################

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




