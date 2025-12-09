# ------------------------------------------------------------
# Utilities
# ------------------------------------------------------------
# Make K distinct random permutations of 1:n (without enumerating n!)
make_support <- function(n, K, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  support <- matrix(NA_integer_, nrow = K, ncol = n)
  seen <- new.env(parent = emptyenv())
  r <- 1L
  while (r <= K) {
    p <- sample.int(n)
    key <- paste(p, collapse = "-")
    if (!exists(key, envir = seen, inherits = FALSE)) {
      support[r, ] <- p
      assign(key, TRUE, envir = seen)
      r <- r + 1L
    }
  }
  support
}

# Calibrate Zipf exponent s so that p1 ≈ target_p1 for given K and q
calibrate_zipf_s <- function(K, target_p1 = 0.4, q = 0, s_lo = 0.1, s_hi = 10) {
  f <- function(s) {
    w <- (1 + q)^(-s)
    denom <- sum(( (1:K) + q )^(-s))
    w / denom - target_p1
  }
  uniroot(f, c(s_lo, s_hi))$root
}

# ------------------------------------------------------------
# Main factory
#   method = "zipf" or "pitmanyor"
#   For Zipf: choose (K, target_p1, q)
#   For Pitman–Yor: choose (d in (0,1), theta > -d), truncated at K
# ------------------------------------------------------------
make_heavy_tail_perm_sampler <- function(
    n,
    K = 100,
    method = c("zipf", "pitmanyor"),
    target_p1 = 0.4,
    zipf_q = 0,
    py_d = 0.5,
    py_theta = 1.0,
    seed = NULL
) {
  method <- match.arg(method)

  # 1) Build K-permutation support
  support <- make_support(n, K, seed = seed)

  # 2) Heavy-tailed probabilities on ranks 1..K
  if (method == "zipf") {
    # Rank permutations by a random shuffle so it's not always the first row
    rank_order <- sample.int(K)
    # Calibrate s so that top mass ~ target_p1
    s <- calibrate_zipf_s(K, target_p1 = target_p1, q = zipf_q)
    w <- ((1:K) + zipf_q)^(-s)
    w <- w / sum(w)
    probs <- numeric(K)
    probs[rank_order] <- w  # assign Zipf weights to a random ranking of support
  } else {
    # Pitman–Yor (GEM(d, theta)) stick-breaking, truncated at K
    # v_k ~ Beta(1 - d, theta + k d), p_k = v_k * prod_{j<k} (1 - v_j)
    v <- numeric(K)
    p <- numeric(K)
    prod_stick <- 1
    for (k in 1:K) {
      a <- 1 - py_d
      b <- py_theta + py_d * (k - 1)
      v[k] <- rbeta(1, a, b)
      p[k] <- prod_stick * v[k]
      prod_stick <- prod_stick * (1 - v[k])
    }
    # (Optionally add the leftover stick to p[K] to ensure exact sum 1)
    p[K] <- p[K] + prod_stick
    # Randomly permute ranks so the top permutation is random in support
    probs <- p[order(runif(K))]
    # If you'd like the expected top mass closer to a target, tweak (py_d, py_theta).
    # Larger d -> heavier tail; smaller theta -> more mass to the head.
  }

  # 3) Sampler closure
  sampler <- function(N = 1) {
    idx <- sample.int(K, size = N, replace = TRUE, prob = probs)
    support[idx, , drop = FALSE]
  }

  list(
    n = n,
    K = K,
    method = method,
    support = support,          # K x n matrix
    probabilities = probs,      # length K, sums to 1
    sample = sampler
  )
}

################################################################################

n <- 5

## ZIPF (power-law) – target top prob around 0.4
S_zipf <- make_heavy_tail_perm_sampler(
  n = n, K = 80, method = "zipf",
  target_p1 = 0.4, zipf_q = 0, seed = 1
)
S_zipf$probabilities          # ~ 0.4 (exact value depends on K & calibration)
which.max(S_zipf$probabilities)    # index of top permutation in S_zipf$support

# Draw samples and estimate empirical probs on the support
samples <- S_zipf$sample(100)
perm_keys <- apply(samples, 1, paste, collapse = "-")
levels_keys <- apply(S_zipf$support, 1, paste, collapse = "-")
emp_zipf <- table(factor(perm_keys, levels = levels_keys)) / nrow(samples)

head(data.frame(
  permutation = levels_keys[order(-S_zipf$probabilities)][1:10],
  true = sort(S_zipf$probabilities, decreasing = TRUE)[1:10],
  empirical = as.numeric(emp_zipf[order(-S_zipf$probabilities)][1:10])
), 10)


################################################################################

permutations <- function(n) {
  if (n == 1) return(matrix(1, nrow = 1))
  prev <- permutations(n - 1)
  out <- NULL
  for (i in 1:n) {
    out <- rbind(out, cbind(i, prev + (prev >= i)))
  }
  colnames(out) <- NULL
  as.matrix(out)
}

if(factorial(n)<500){
  centering_grid <- permutations(n)
}else{
  centering_grid <- t(replicate(500, sample.int(n)))
}
#centering_grid <- t(replicate(500, sample.int(n)))  # grid = centering set
#centering_grid <- rbind(c(1,2,3), c(1,3,2), c(2,3,1),c(2,1,3), c(3,1,2), c(3,2,1))

kef_result <- kef(samples = samples,grids = centering_grid,lambda = 1,tau = 1/100000,
    data_type = "order")



## ------------------------------------------------------------
## Map everything to each sampled permutation (row-wise)
## ------------------------------------------------------------

# character key for each row of the *support*
levels_keys <- apply(S_zipf$support, 1, paste, collapse = "-")

# character key for each sampled permutation
perm_keys <- apply(samples, 1, paste, collapse = "-")

# empirical probabilities on the support (you already computed this)
emp_zipf <- table(factor(perm_keys, levels = levels_keys)) / nrow(samples)

# indices of each sampled permutation in the support
idx_support <- match(perm_keys, levels_keys)

# build data frame:
kef_results_df <- data.frame(
  sample              = perm_keys,                       # permutation as string
  true_prob           = S_zipf$probabilities[idx_support],
  empirical_prob      = as.numeric(emp_zipf[idx_support]),
  kef_estim_dens      = kef_result$dens_sample          # kef estimated prob (per sample)
)

# build data frame:
kef_results_df_grid <- data.frame(
  grid                = apply(centering_grid, 1, paste, collapse = "-"),                     # permutation as string
  kef_estim_dens      = kef_result$dens_grids         # kef estimated prob (per sample)
)



#sum(get_base_measures(samples,data_type = "order"))
#get_base_measures(samples,data_type = "order")

head(kef_results_df,10)

library(dplyr)
kef_results_df_unique<- kef_results_df %>% unique()
kef_results_df_unique$kef_normalized_prob <- 1/sum(kef_results_df_unique$kef_estim_dens) * kef_results_df_unique$kef_estim_dens
kef_results_df_unique

kef_results_df <- merge(kef_results_df,kef_results_df_unique,by = names(kef_results_df))
#sum(kef_results_df$empirical_prob)
#sum(kef_results_df$kef_estim_prob)
View(kef_results_df)

base_g <- get_base_measures(centering_grid, "order")
base_g
base_s <-  get_base_measures(samples, "order")
base_s
