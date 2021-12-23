parameter_combiner <- function(method, data, adjacency_matrix, graph_configs = NULL)
{
  # unpack necessary parameters
  method_type <- method$method_type
  hyperparameter_configs <- method$hyperparameter_configs
  
  if(method_type == 'appr')
  {
    parameter_combinations <- appr_parameter_combiner(hyperparameter_configs, data, adjacency_matrix)
  } else if(method_type == 'local_density')
  {
    parameter_combinations <- local_density_parameter_combiner(hyperparameter_configs, data, adjacency_matrix,
                                                         graph_configs)
  } else if(method_type == 'spectral')
  {
    parameter_combinations <- spectral_parameter_combiner(hyperparameter_configs, data)
  } else if(method_type == 'ppr')
  {
    parameter_combinations <- ppr_parameter_combiner(hyperparameter_configs, data)
  } else if(method_type == 'conductance')
  {
   parameter_combinations <- conductance_parameter_combiner(hyperparameter_configs,adjacency_matrix,
                                                            data)
  } else if(method_type == 'density')
  {
    parameter_combinations <- density_parameter_combiner(hyperparameter_configs, data)
  }
  return(parameter_combinations)
}

density_parameter_combiner <- function(density_hyperparameter_configs, data)
{
  # Input: density_hyperparameter_configs
  #        - $density_cut (list)
  #        - - $type ('tune')
  
  # unpack necessary parameters
  density_cut_type <- density_hyperparameter_configs$density_cut$type
  density_cut_length <- density_hyperparameter_configs$density_cut$length
  N <- density_hyperparameter_configs$N
  n <- nrow(data)
  
  if(density_cut_type == 'tune')
  {
    density_cut <- seq(1, N + n, length.out = density_cut_length)
  }
  n_parameters <- 2 # density_cut, N extra points
  parameter_combinations <- matrix(ncol = n_parameters, nrow = length(density_cut),
                                   dimnames = list(NULL, c('density_cut', 'N')))
  for(ii in 1:length(density_cut))
  {
    parameter_combinations[ii,] = c(density_cut[ii], N)
  }
  return(parameter_combinations)
}

appr_parameter_combiner <- function(appr_hyperparameter_configs, data, adjacency_matrix)
{
  # Input: appr_hyperparameter_configs (list of parameters)
  #        - $seed_node_location (vector) (location of seed node)
  #        - $teleportation_parameter (list)
  #        - - $type (flag) ('tune')
  #        - - $range (2-vector) (upper and lower bounds for teleportation vector)
  #        - - $length (numeric) number of different teleportation parameters we will try
  #        - $approximation_parameter (list) 
  #        - - $type (flag) ('tune)
  #        - - $proportion_volume (numeric) volume of the largest cluster algorithm might recover
  
  # unpack parameters
  seed_node_manifold_location <- appr_hyperparameter_configs$seed_node_manifold_location
  teleportation_parameter_range <- appr_hyperparameter_configs$teleportation_parameter$range
  teleportation_parameter_length <- appr_hyperparameter_configs$teleportation_parameter$length
  approximation_parameter_proportion_volume <- appr_hyperparameter_configs$approximation_error$proportion_volume
  
  
  # calculate grids for tuning parameters
  d <- ncol(data)
  seed_node_location <- intrinsic_to_extrinsic(d,seed_node_manifold_location)
  v_grid <- which.min(apply(data, 1, FUN = function(x){sum((x - seed_node_location)^2)}))
  
  if(appr_hyperparameter_configs$teleportation_parameter$type == 'tune')
  {
    alpha_grid <- exp( seq(log(teleportation_parameter_range[1]), 
                           log(teleportation_parameter_range[2]), 
                           length.out = teleportation_parameter_length ) )
  }
  
  if(appr_hyperparameter_configs$approximation_parameter$type == 'tune')
  {
    vol_G <- sum(adjacency_matrix)
    vol_0_grid <- 2^((1:ceiling(log2(vol_G/4))))
    eps_grid <- (1 / 10) * (1 / vol_0_grid) # choice of 1/10 motivated by theory (Zhu 2013)
  }
  n_parameters <- 3 # v, alpha, eps
  n_combinations <- length(v_grid) * length(alpha_grid) * length(eps_grid)
  parameter_combinations <- matrix(ncol = n_parameters, nrow = n_combinations)
  colnames(parameter_combinations) <- c('seed_node','teleportation_parameter', 'approximation_parameter')
  ii = 1 # row index
  for(v in v_grid)
  {
    for(alpha in alpha_grid)
    {
      for(eps in eps_grid)
      {
        parameter_combinations[ii,] <- c(v,alpha,eps)
        ii = ii + 1 
      }
    }
  }
  return(parameter_combinations)
}

local_density_parameter_combiner <- function(density_hyperparameter_configs, data, adjacency_matrix,
                                       graph_configs)
{
  # unpack necessary parameters
  seed_node_manifold_location <- density_hyperparameter_configs$seed_node_manifold_location
  graph_type <- graph_configs$graph_type
  connectivity <- graph_configs$connectivity

  # calculate grids for tuning parameters
  d <- ncol(data)
  seed_node_location <- intrinsic_to_extrinsic(d,seed_node_manifold_location)
  v_grid <- which.min(apply(data, 1, FUN = function(x){sum((x - seed_node_location)^2)}))
  if(graph_type == 'kNN')
  {
    prune_by <- 'kth_neighbor_radius'
  } else if (graph_type == 'epsilon')
  {
    prune_by <- 'degree'
  }
  parameter_combinations <- matrix(c(v_grid, connectivity, prune_by), 
                                   ncol = 3,
                                   dimnames = list(NULL, c('seed_node', 'connectivity','prune_by')))
  return(parameter_combinations)
}

conductance_parameter_combiner <- function(conductance_hyperparameter_configs, adjacency_matrix, data)
{
  # unpack necessary parameters
  seed_node_manifold_location <- conductance_hyperparameter_configs$seed_node_manifold_location
  tol <- conductance_hyperparameter_configs$tol
  c <- conductance_hyperparameter_configs$c
  f0 <- conductance_hyperparameter_configs$f0
  
  # calculate grids for tuning parameters
  d <- ncol(data)
  seed_node_location <- intrinsic_to_extrinsic(d,seed_node_manifold_location)
  v_grid <- which.min(apply(data, 1, FUN = function(x){sum((x - seed_node_location)^2)}))
  
  parameter_combinations <- matrix(c(v_grid,tol,c,f0), ncol = 4, dimnames = list(NULL, 
                                                                                 c('seed_node','tol','c','f0')))
  return(parameter_combinations)
}

spectral_parameter_combiner <- function(spectral_hyperparameter_configs, data)
{
  # unpack necessary parameters
  seed_node_manifold_location <- spectral_hyperparameter_configs$seed_node_manifold_location
  
  # calculate grids for tuning parameters
  d <- ncol(data)
  seed_node_location <- intrinsic_to_extrinsic(d,seed_node_manifold_location)
  v_grid <- which.min(apply(data, 1, FUN = function(x){sum((x - seed_node_location)^2)}))
  
  parameter_combinations <- matrix(v_grid, ncol = 1, dimnames = list(NULL, c('seed_node')))
  return(parameter_combinations)
}

ppr_parameter_combiner <- function(ppr_hyperparameter_configs, data)
{
  # Input: ppr_hyperparameter_configs (list of parameters)
  #        - $seed_node_manifold_location (vector) location of seed node on manifold
  #        - $teleportation_parameter (list)
  #        - - $type (flag) ('tune')
  #        - - $range (2-vector) (upper and lower bounds for teleportation vector)
  #        - - $length (numeric) number of different teleportation parameters we will try
  
  # unpack parameters
  seed_node_manifold_location <- ppr_hyperparameter_configs$seed_node_manifold_location
  teleportation_parameter_range <- ppr_hyperparameter_configs$teleportation_parameter$range
  teleportation_parameter_length <- ppr_hyperparameter_configs$teleportation_parameter$length
  
  
  # calculate grids for tuning parameters
  d <- ncol(data)
  seed_node_location <- intrinsic_to_extrinsic(d,seed_node_manifold_location)
  v_grid <- which.min(apply(data, 1, FUN = function(x){sum((x - seed_node_location)^2)}))
  
  if(ppr_hyperparameter_configs$teleportation_parameter$type == 'tune')
  {
    alpha_grid <- exp( seq(log(teleportation_parameter_range[1]), 
                           log(teleportation_parameter_range[2]), 
                           length.out = teleportation_parameter_length ) )
  }
  
  n_parameters <- 2
  parameter_combinations <- matrix(nrow = length(v_grid) * length(alpha_grid), 
                                   ncol = n_parameters, 
                                   dimnames = list(NULL, c('seed_node','teleportation_parameter')))
  ii = 1
  for(v in v_grid)
  {
    for(alpha in alpha_grid)
    {
      parameter_combinations[ii,] <- c(v, alpha)
      ii = ii + 1
    }
  }
  return(parameter_combinations)
}