# sample size
n <- 800

# distribution
P <- list(distribution_class = 'two_moons',
          d = 4,
          centers = matrix(c(-.5,-.075,
                             0,.075),ncol = 2, byrow = T),
          r = .5,
          sigma = .075,
          cluster_probabilities = c(.5,.5),
          theta_dist = 'uniform'
) 

# graph
graph_configs <- list(graph_type = 'kNN',
                      connectivity = 6)

# clustering methods
cluster_methods <- list(
  method_D = list(
    method_type = 'conductance',
    hyperparameter_configs = list(
      seed_node_manifold_location = c(.05, -.35),
      tol = 10^(-6),
      c = 1,
      f0 = 'tune_spectral'
    )
  ),
  method_A = list(
    method_type = 'ppr',
    hyperparameter_configs = list(
      seed_node_manifold_location = c(.05, -.425),
      teleportation_parameter = list(
        type = 'tune',
        range = c(.001,.3,25),
        length = 25
      )
    )
  ),
  method_B = list(
    method_type = 'local_density',
    hyperparameter_configs = list(
      seed_node_manifold_location = c(.05, -.35)
    )
  ),
  method_C = list(
    method_type = 'spectral',
    hyperparameter_configs = list(
      seed_node_manifold_location = c(.05, -.35)
    )
  )
)

evaluation_configs <- list(metric = 'ncut')


# cluster_methods <- list(
#   method_A = list(
#     graph_type = 'kernel',
#     kernel_type = 'uniform',
#     h = 'quantile',
#     quantile = .2,
#     cluster_type = 'PPR',
#     seed_node = 'nearest_origin',
#     teleportation_parameter = list(
#       type = 'tune',
#       L = .001,
#       U = .3,
#       length = 25),
#     target_volume = 'tune'
#   ),
#   method_B = list(
#     graph_type = 'kernel',
#     kernel_type = 'uniform',
#     h = 'quantile',
#     quantile = .1,
#     cluster_type = 'PPR',
#     seed_node = 'nearest_origin',
#     teleportation_parameter = list(
#       type = 'tune',
#       L = .001,
#       U = .3,
#       length = 25),
#     target_volume = 'tune'
#   ),
#   method_C = list(
#     graph_type = 'kernel',
#     kernel_type = 'uniform',
#     h = 'quantile',
#     quantile = .05,
#     cluster_type = 'PPR',
#     seed_node = 'nearest_origin',
#     teleportation_parameter = list(
#       type = 'tune',
#       L = .001,
#       U = .3,
#       length = 25),
#     target_volume = 'tune'
#   ),
#   method_D = list(
#     graph_type = 'kernel',
#     kernel_type = 'uniform',
#     h = 'quantile',
#     quantile = .01,
#     cluster_type = 'PPR',
#     seed_node = 'nearest_origin',
#     teleportation_parameter = list(
#       type = 'tune',
#       L = .001,
#       U = .3,
#       length = 25),
#     target_volume = 'tune'
#   )
# )

# simulation hyperparameters (?)
n_sims <- 100