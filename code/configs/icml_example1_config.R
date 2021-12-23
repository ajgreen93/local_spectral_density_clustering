# sample size
n <- 5000

# distribution
P <- list(distribution_class = 'gaussian_mixtures',
          d = 1,
          cluster_probabilities = c(.5,.5),
          mu_type = 'evenly_spaced',
          Sigma_type = 'scaled',
          Sigma_scale = 1/2)

# density cluster parameters
density_cluster_parameters <- list(
  lambda_path = seq(1.5,1, length.out = 5),
  sigma = .01,
  r = 'sigma/4d'
)
