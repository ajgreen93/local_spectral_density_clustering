local_clusterer <- function(data, adjacency_matrix, method, parameters)
{
  # Input: similarity (n vector)
  #        adjacency_matrix (n x n matrix)
  #        evaluation_configs (list)
  #        -$metric (flag) ('ncut')
  
  # unpack parameters
  method_type <- method$method_type
  evaluation_configs <- method$evaluation_configs
  
  # Compute similarities
  if(method_type %in% c('ppr', 'local_density', 'conductance', 'spectral'))
  {
    similarity <- similarity_scorer(method_type, adjacency_matrix, parameters, data)
    cluster <- sweep_cut(similarity, evaluation_configs, adjacency_matrix)
  } 
  return(list(similarity = similarity,
              cluster = cluster))
}

tree_clusterer <- function(data, adjacency_matrix, method, parameters)
{
  # Input: data (n x d matrix)
  #        adjacency_matrix (n x n matrix)
  #        method (list of configurations which define a method)
  #        parameters (inputs to method)
  #        N (number samples to augment with)
  
  # unpack necessary parameters
  method_type <- method$method_type
  
  if(method_type == 'density')
  {
    cluster_tree <- density_cluster_tree(data, adjacency_matrix, parameters)$empirical_cluster_tree
  }
  
}

tree_to_clustering <- function(cluster_tree, parameters, P)
{
  # unpack necessary parameters
  p <- length(P$cluster_probabilities)
  
  # cut the (density) cluster tree 
  clustering <- cluster_tree[parameters['density_cut'],]
  
  # all nodes with density less than the cut should not be in a cluster
  cluster_sizes <- sort(table(clustering), decreasing = T)
  cluster_names <- as.numeric(names(cluster_sizes)[cluster_sizes  > 1])
  clustering[!(clustering %in% cluster_names)] <- 0
  
  return(clustering)
}
density_cluster_tree <- function(data,adjacency_matrix, parameters)
{
  density_data <- density(NULL, data, 0) # TODO: make P optional input to density
  cluster_tree_data <- cluster_tree(NULL, data, approx = F, density_data)
  return(cluster_tree_data)
}

sweep_cut <- function(similarity, evaluation_configs, adjacency_matrix)
{
  # Initialize best case cluster and evaluation
  cluster0 <- c()
  evaluation0 <- 1
  
  # Compute best cluster
  similarity_ranking <- order(similarity,decreasing = T)
  for(jj in 1:length(similarity_ranking))
  {
    cluster <- similarity_ranking[1:jj]
    evaluation <- evaluator(evaluation_configs,cluster,adjacency_matrix)$evaluation
    
    # keep track of the best cluster
    if(evaluation < evaluation0)
    {
      cluster0 <- cluster
      evaluation0 <- evaluation
    }
  }
  
  return(cluster0)
}

ppr_cluster <- function(A, X, cluster_method)
{
  # choose optimal conductance cluster over possible options
  target_set_conductance <- 1
  cluster <- c()
  for(v in v_grid)
  {
    for(alpha in alpha_grid)
    {
      for(vol_0 in vol_0_grid)
      {
        
        p <- appr(A, v, alpha, eps)
        C <- best_sweep_cut(A, p, vol_0)
        phi <- conductance(A, C)
        if(phi < target_set_conductance)
        {
          cluster <- C
          target_set_conductance <- phi
          print(target_set_conductance)
        }
        # print(paste('v = ', v, 'alpha = ', alpha, 'vol_0 = ', vol_0), collapse = T)
      }
    }
  }
  return(cluster)
}

#------------ OBSOLETE ----------------#
best_sweep_cut <- function(A, p, vol_0)
{
  # set up parameters
  phi <- 1
  cluster <- c()
  
  # compute normalized PPR vector
  degree <- rowSums(A)
  q <- p/degree
  ordering <- order(q, decreasing = TRUE)
  q_sort <- q[ordering]
  
  # find best sweep cut according to conductance, with boundary conditions on search due to (Zhu)
  candidates <- which(q_sort >= 1/5 * 1 / vol_0 & q_sort <= 3/5 * 1 / vol_0)
  if(length(candidates) == 0)
  {
    cluster = c()
  } else{
    for(candidate in candidates)
    {
      C_candidate <- ordering[1:candidate]
      phi_candidate <- conductance(A, C_candidate)
      if(phi_candidate < phi)
      {
        phi <- phi_candidate
        cluster <- C_candidate
      }
    }
  }
  return(cluster)
}

