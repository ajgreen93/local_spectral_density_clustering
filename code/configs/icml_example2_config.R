# sample size
n <- 400

# distribution
P <- list(distribution_class = 'gaussian_mixtures',
          d = 1,
          cluster_probabilities = c(.5,.5),
          mu_type = 'evenly_spaced',
          Sigma_type = 'scaled',
          Sigma_scale = 1) 

# clustering methods
cluster_methods <- list(
  method_A = list(
    graph_type = 'kernel',
    kernel_type = 'uniform',
    h = 'quantile',
    quantile = .1,
    cluster_type = 'PPR',
    seed_node = 'nearest_origin',
    teleportation_parameter = list(
      type = 'tune',
      L = .001,
      U = .3,
      length = 25),
    target_volume = 'tune'
  )
)

# simulation hyperparameters (?)
n_sims <- 100