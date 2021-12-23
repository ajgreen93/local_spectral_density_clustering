# distribution
P <- list(distribution_class = 'gaussian_mixtures',
                   d = 2,
                   cluster_probabilities = c(.5,.5),
                   mu_type = 'evenly_spaced',
                   Sigma = matrix(c(33/32,-31/32,
                                    -31/32,33/32),nrow = 2)/4,
                   Sigma_type = 'None') 

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

density_cluster_parameters <- list(
  lambda_path = seq(1.5,1, length.out = 5),
  sigma = .12*4,
  r = 'sigma/4d'
)

# simulation hyperparameters (?)
n_sims <- 5