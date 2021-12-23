clusterer <- function(X, cluster_methods)
{
  # objects for storing output of cluster methods, such as local clusters
  cluster_data <- vector(mode = 'list',length(cluster_methods))
  names(cluster_data) = names(cluster_methods)
  
  # run each cluster method
  for(cluster_method_name in names(cluster_methods) )
  {
    cluster_method <- cluster_methods[[cluster_method_name]]
    
    # form graph
    graph_type <- cluster_method$graph_type
    if(graph_type == 'kNN')
    {
      A <- knn_grapher(X, cluster_method)
    } else if(graph_type == 'kernel')
    {
      A <- kernel_grapher(X, cluster_method)
    }
    
    # calculate cluster
    cluster_type <- cluster_method$cluster_type
    if(cluster_type == 'PPR')
    {
      cluster <- ppr_cluster(A, X, cluster_method)
    }
    cluster_data[[cluster_method_name]] <- list(
      cluster_method = cluster_method,
      A = A,
      cluster = cluster
    )
  }
  return(cluster_data)
}





graphical_evaluator_example2_plot1 <- function(X, cluster_data)
{
  component_densities <- gaussian_mixture_densities(P)
  f_nk <- matrix(nrow = nrow(X), ncol = length(component_densities))
  for(i in 1:nrow(X))
  {
    f_nk[i,] <- sapply(component_densities,FUN = function(density){density(X[i,])})
  }
  density_cluster_membership <- apply(f_nk,1,FUN = function(f){which.max(f)})
  cluster_probabilities <- P$cluster_probabilities
  f_n <- f_nk %*% cluster_probabilities
  
  origin_node <- which.min(apply(X,1,FUN = function(x){sum(x^2)}))
  density_cluster <- which(density_cluster_membership == density_cluster_membership[origin_node])
  
  cluster_matrix <- matrix(0,nrow = nrow(X), ncol = length(cluster_data))
  colnames(cluster_matrix) <- names(cluster_data)
  for(cluster_method_name in names(cluster_data))
  {
    cluster_results <- cluster_data[[cluster_method_name]]
    estimated_cluster <- cluster_results$cluster
    cluster_matrix[estimated_cluster,cluster_method_name] <- 1
  }
  
  plot_df <- as.data.frame(cluster_matrix)
  plot_df['X'] <- X[,1]
  plot_df['Y'] <- f_n
  plot_df['cluster_membership'] <- density_cluster_membership
  
  example2_plot1 <- ggplot(data = plot_df, aes(x = X, y = Y)) + 
    geom_point(aes(colour = as.factor(method_A), shape = as.factor(cluster_membership))) + 
    scale_shape_manual(name = 'Mixture component',
                       values = c(2,3,1), 
                       labels = c('Cluster 1', 'Cluster 2', 'Cluster 3')) +
    scale_color_manual(name = 'PPR cluster',values = c('red3', 'blue3'),
                       labels = c('Out of PPR cluster', 'In PPR cluster')) +
    xlab('x') + ylab('density') +
    theme_classic()
}

graphical_evaluator_example2_plot2 <- function(P, X, cluster_data)
{
  component_densities <- gaussian_mixture_densities(P)
  f_nk <- matrix(nrow = nrow(X), ncol = length(component_densities))
  for(i in 1:nrow(X))
  {
    f_nk[i,] <- sapply(component_densities,FUN = function(density){density(X[i,])})
  }
  density_cluster_membership <- apply(f_nk,1,FUN = function(f){which.max(f)})
  
  origin_node <- which.min(apply(X,1,FUN = function(x){sum(x^2)}))
  density_cluster <- which(density_cluster_membership == density_cluster_membership[origin_node])
  
  cluster_matrix <- matrix(0,nrow = nrow(X), ncol = length(cluster_data))
  colnames(cluster_matrix) <- names(cluster_data)
  for(cluster_method_name in names(cluster_data))
  {
    cluster_results <- cluster_data[[cluster_method_name]]
    estimated_cluster <- cluster_results$cluster
    cluster_matrix[estimated_cluster,cluster_method_name] <- 1
  }
  
  
  plot_df <- as.data.frame(cluster_matrix)
  plot_df['X'] <- X[,1]
  plot_df['Y'] <- X[,2]
  plot_df['cluster_membership'] <- density_cluster_membership
  
  example2_plot2 <- ggplot(data = plot_df, aes(x = X, y = Y)) + 
    geom_point(aes(colour = as.factor(method_A), shape = as.factor(density_cluster_membership))) + 
    scale_shape_manual(name = 'Mixture component',
                       values = c(2,3), 
                       labels = c('Cluster 1', 'Cluster 2')) +
    scale_color_manual(name = 'PPR cluster',values = c('red3', 'blue3'),
                       labels = c('Out of PPR cluster', 'In PPR cluster')) +
    xlab('x1') + ylab('x2') + coord_fixed(xlim = c(-.5,.75), ylim = c(-.5,.75)) + 
    theme_classic()
}

graphical_evaluator_plot3 <- function(P_sequence,metric_evaluation_table)
{
  average_symmetric_set_difference_performance <- rowMeans(metric_evaluation_table)
  cluster_center_distance <- P_sequence$mu_distance_scale * .5
  
  plot_df = data.frame(X = cluster_center_distance, Y = average_symmetric_set_difference_performance)
  plot3 <- ggplot(data = plot_df, aes(x = X, y = Y)) + geom_line(color = 'red') +
    geom_point(color = 'red', shape = 2) + 
    xlab('Distance between component centers') + ylab('symmetric set difference') + 
    theme_classic()
}

graphical_evaluator_example2_plot3<- function(P, X, cluster_data)
{
  component_densities <- gaussian_mixture_densities(P)
  f_nk <- matrix(nrow = nrow(X), ncol = length(component_densities))
  for(i in 1:nrow(X))
  {
    f_nk[i,] <- sapply(component_densities,FUN = function(density){density(X[i,])})
  }
  density_cluster_membership <- apply(f_nk,1,FUN = function(f){which.max(f)})
  
  origin_node <- which.min(apply(X,1,FUN = function(x){sum(x^2)}))
  density_cluster <- which(density_cluster_membership == density_cluster_membership[origin_node])
  
  cluster_matrix <- matrix(0,nrow = nrow(X), ncol = length(cluster_data))
  colnames(cluster_matrix) <- names(cluster_data)
  for(cluster_method_name in names(cluster_data))
  {
    cluster_results <- cluster_data[[cluster_method_name]]
    estimated_cluster <- cluster_results$cluster
    cluster_matrix[estimated_cluster,cluster_method_name] <- 1
  }
  
  
  plot_df <- as.data.frame(cluster_matrix)
  plot_df['X'] <- X[,1]
  plot_df['Y'] <- X[,2]
  plot_df['cluster_membership'] <- density_cluster_membership
  
  example2_plot3 <- ggplot(data = plot_df, aes(x = X, y = Y)) + 
    geom_point(aes(colour = as.factor(method_A), shape = as.factor(density_cluster_membership))) + 
    scale_shape_manual(name = 'Mixture component',
                       values = c(2,3), 
                       labels = c('Cluster 1', 'Cluster 2')) +
    scale_color_manual(name = 'PPR cluster',values = c('red3', 'blue3'),
                       labels = c('Out of PPR cluster', 'In PPR cluster')) +
    xlab('x1') + ylab('x2') + coord_fixed(xlim = c(-1.5,1.5), ylim = c(-1.5,1.5)) + 
    theme_classic()
}

graphical_evaluator_example2_plot4 <- function(P,X,cluster_data)
{
  component_densities <- gaussian_mixture_densities(P)
  f_nk <- matrix(nrow = nrow(X), ncol = length(component_densities))
  for(i in 1:nrow(X))
  {
    f_nk[i,] <- sapply(component_densities,FUN = function(density){density(X[i,])})
  }
  density_cluster_membership <- apply(f_nk,1,FUN = function(f){which.max(f)})
  cluster_probabilities <- P$cluster_probabilities
  f_n <- f_nk %*% cluster_probabilities
  
  origin_node <- which.min(apply(X,1,FUN = function(x){sum(x^2)}))
  density_cluster <- which(density_cluster_membership == density_cluster_membership[origin_node])
  
  cluster_matrix <- matrix(0,nrow = nrow(X), ncol = length(cluster_data))
  colnames(cluster_matrix) <- names(cluster_data)
  for(cluster_method_name in names(cluster_data))
  {
    cluster_results <- cluster_data[[cluster_method_name]]
    estimated_cluster <- cluster_results$cluster
    cluster_matrix[estimated_cluster,cluster_method_name] <- 1
  }
  
  plot_df <- as.data.frame(cluster_matrix)
  plot_df['X'] <- X[,1]
  plot_df['Y'] <- f_n
  plot_df['cluster_membership'] <- density_cluster_membership
  
  example2_plot4 <- ggplot(data = plot_df, aes(x = X, y = Y)) + 
    geom_point(aes(colour = as.factor(method_A), shape = as.factor(cluster_membership))) + 
    scale_shape_manual(name = 'Mixture component',
                       values = c(2,3,1), 
                       labels = c('Cluster 1', 'Cluster 2', 'Cluster 3')) +
    scale_color_manual(name = 'PPR cluster',values = c('red3', 'blue3'),
                       labels = c('Out of PPR cluster', 'In PPR cluster')) +
    xlab('x') + ylab('density') +
    theme_classic()
}
metric_evaluator <- function(P,X,cluster_data,density_cluster_parameters)
{
  component_densities <- gaussian_mixture_densities(P)
  f_nk <- matrix(nrow = nrow(X), ncol = length(component_densities))
  for(i in 1:nrow(X))
  {
    f_nk[i,] <- sapply(component_densities,FUN = function(density){density(X[i,])})
  }
  density_cluster_membership <- apply(f_nk,1,FUN = function(f){which.max(f)})
  
  origin_node <- which.min(apply(X,1,FUN = function(x){sum(x^2)}))
  density_cluster <- which(density_cluster_membership == density_cluster_membership[origin_node])
  
  hatC <- cluster_data$method_A$cluster
  symmetric_set_difference <- max(length(setdiff(density_cluster,hatC))/length(density_cluster),
                                  length(setdiff(hatC,density_cluster))/length(hatC))
  
  # for a range of lambda
  # lambda_seq = seq(0,3,length.out = 150)
  # symmetric_set_difference = 1
  # for(lambda in lambda_seq)
  # {
  #   # take the local density cluster
  #   C_lambda <- density_cluster(P,X,density_cluster_parameters, lambda)
  #   hatC <- cluster_data$method_A$cluster
  #   candidate_symmetric_set_difference = max(length(setdiff(C_lambda,hatC))/length(C_lambda),
  #                                            length(setdiff(hatC,C_lambda))/length(hatC))
  #   # print(candidate_symmetric_set_difference)
  #   if(candidate_symmetric_set_difference < symmetric_set_difference)
  #   {
  #     symmetric_set_difference <- candidate_symmetric_set_difference
  #   }
  # }
  return(symmetric_set_difference)
}