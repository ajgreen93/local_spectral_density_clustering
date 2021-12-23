evaluator <- function(evaluation_configs, cluster, adjacency_matrix, data, P)
{
  # unpack necessary parameters
  metric <- evaluation_configs$metric
  
  # compute the desired metric
  if(metric == 'ncut')
  {
    evaluation_data <- list(
      evaluation = normalized_cut(cluster,adjacency_matrix)
    )
  } else if(metric == 'density_clustertree')
  {
    evaluation_data <- density_clustertree_recovery(cluster, data, P)
  } else if(metric == 'oracle_density_clustering')
  {
    evaluation_data <- oracle_density_evaluation(cluster, P)
  }
  
  return(evaluation_data)
}

oracle_density_evaluation <- function(clustering, P)
{
  # oracle_density_evaluation should give lower (better) scores
  # to clusterings which assign most points to the top p cluster
  # where p is the number of mixture components in P
  # 
  # INPUT: clustering (n x 1 vector of cluster membership for each data points)
  #        P (list of specifications for mixture model)
  
  # unpack necessary parameters
  p <- length(P$cluster_probabilities)
  
  # size of pth largest cluster
  cluster_sizes <- sort(table(clustering[clustering != 0]), decreasing = T)
  if(length(cluster_sizes) < p)
  {
    pth_largest_cluster_size <- 1
  } else{
    pth_largest_cluster_size <- cluster_sizes[p]
  }
  
  return(list(
    evaluation = 1/pth_largest_cluster_size
    )
  )
  
}

density_clustertree_recovery <- function(cluster, data, P, cluster_tree_data,
                                         evaluation_metric = 'set_difference')
{
  # TODO: function definition here
  # Output: (list)
  #         - C_nlprime: indices of data points with density at least
  #                      that of threshold for density cluster but not 
  #                      part of best density cluster
  
  if(is.null(cluster_tree_data))
  {
    cluster_tree_data <- cluster_tree(P, data)
  }
  empirical_cluster_tree <- cluster_tree_data$empirical_cluster_tree
  density_rank <- cluster_tree_data$density_rank
  n_tree_cuts <- nrow(empirical_cluster_tree)
  
  # For each radius value, find connected components in the graph
  evaluation0_min <- 1
  evaluation0 <- rep(1,n_tree_cuts)
  for(ii in 1:n_tree_cuts )
  {
    U <- empirical_cluster_tree[ii,]
    U[table(U)[U] == 1] = 0
    cluster_labels <- setdiff(unique(U),c(0))
    for(l in cluster_labels)
    {
      C_nl <- which(U == l)
      C_nlprime <- which(U != l & U > 0)
      if(evaluation_metric == 'set_difference')
      {
        evaluation <- length(setdiff(C_nl,cluster)) / length(C_nl) +
                      length(setdiff(cluster,C_nl)) / length(cluster)
      } else{
        evaluation <- max( length(setdiff(C_nl,cluster)) / length(C_nl),
                           ifelse(length(C_nlprime) > 0,
                                  1 - length(setdiff(C_nlprime,cluster)) / length(C_nlprime),
                                  0) )
      }
      
      # when two density clusters are equally well recovered, default to the larger one
      if(evaluation <= evaluation0[ii])
      {
        evaluation0[ii] <- evaluation
        if(evaluation <= evaluation0_min)
        {
          evaluation0_min <- evaluation
          C_n0 <- C_nl
          C_n0prime <- C_nlprime
        }
      }
    }
    print(ii)
  }
  return(
    list(evaluation = evaluation0_min,
         cluster_tree_evaluation = matrix(c(evaluation0,density_rank), 
                                          ncol = 2,
                                          dimnames = list(NULL,c('cluster_recovery','density_threshold'))),
         closest_density_cluster = C_n0,
         other_density_clusters = C_n0prime)
  )
}