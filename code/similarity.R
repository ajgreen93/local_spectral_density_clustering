similarity_scorer <- function(method_type, adjacency_matrix, parameters, data)
{
  # Input: method_type (flag) local clustering method ('appr', 'density', 'spectral', 'ppr')
  #        adjacency_matrix (n x n matrix)
  #        parameters (named vector) vector of parameters for use
  
  if(method_type == 'appr')
  {
    similarity <- appr(adjacency_matrix,parameters)
  } else if(method_type == 'local_density')
  {
    similarity <- local_density(adjacency_matrix,parameters,data)
  } else if(method_type == 'spectral')
  {
    similarity <- spectral(adjacency_matrix,parameters)
  } else if(method_type == 'ppr')
  {
    similarity <- ppr(adjacency_matrix,parameters)
  } else if(method_type == 'conductance')
  {
    similarity <- conductance(adjacency_matrix, parameters)
  }
  return(similarity)
}

appr <- function(A, parameters)
{
  # Input: A (n x n matrix) adjacency matrix
  #        parameters (vector)
  #        -['seed_node']
  #        -['teleportation_parameter']
  #        -['approximation_parameter']
  
  # unpack necessary parameters
  v <- as.integer(parameters['seed_node']) # integer conversion necessary for Julia call
  alpha <- parameters['teleportation_parameter']
  eps <- parameters['approximation_parameter']
  
  # compute personalized PageRank
  A_list <- which(A == 1, arr.ind = T)
  aprJL <- julia_call("apr_wrapper", A_list, v, alpha, eps)
  appr_vector <- unlist(aprJL) * rowSums(A)[as.numeric(names(unlist(aprJL)))]
  ppr_vector <- numeric(nrow(A)); ppr_vector[as.numeric(names(appr_vector))] <- appr_vector
  
  similarity <- ppr_vector / rowSums(A)
  
  return(similarity)
}

conductance <- function(A, parameters)
{
  # Input: A (n x n matrix) adjacency matrix,
  #        data (n x d matrix) data matrix
  # Output: second eigenvector as 'similarity' vector
  
  # unpack necessary parameters
  seed_node <- strtoi(parameters['seed_node'])
  tol <- as.numeric(parameters['tol'])
  f0 <- parameters['f0']
  c <- strtoi(parameters['c'])
  
  if(f0 == 'tune_spectral')
  {
    f0 <- similarity_scorer('spectral', adjacency_matrix, parameters = c('seed_node' = seed_node) )
  }
  
  n <- nrow(A)
  
  f_new <- f0 - median(f0)
  
  plot(data[f0 > 0,])
  energy_old <- energy(A,f0) + tol + 1
  energy_new <- energy(A,f_new)
  while(energy_old - energy_new > tol)
  {
    ## 0. new becomes old
    energy_old <- energy_new
    f <- f_new
    
    ## 1. identify element of subgradient
    v <- numeric(length(f))
    v[f > 0] = 1
    v[f < 0] = -1
    v[f == 0] = (sum(f < 0) - sum(f > 0)) / (sum(f == 0))
    
    ## 2. move f in the direction of the subgradient
    g <- f + c * v
    
    ## 3. minimize total variation in the area around g
    D <- incidence(A)
    m <- nrow(D)
    model <- list()
    model$modelsense <- 'min'
    # TODO: make this faster by eliminating the use of rbind
    model$Q <- energy(A, f) / (2*c) * rbind(
      cbind(
        matrix(0, ncol = m, nrow = m),
        matrix(0, ncol = n, nrow = m)
      ),
      cbind(
        matrix(0, ncol = m, nrow = n),
        diag(n)
      )
    )
    model$obj <- c(rep(1,m), - energy(A, f) / (c) * g)
    model$A <- rbind(
      cbind(
        - diag(m),
        D
      ),
      cbind(
        -diag(m),
        - D
      )
    )
    model$rhs <- rep(0, 2*m)
    
    result <- gurobi(model)
    h <- result$x[(m + 1):(m + n)]
    
    # 4. center
    h <- h - median(h)
    
    # 5. scale
    f_new <- h / norm(h,0)
    energy_new <- energy(A,f_new)
    print(energy_new)
    plot(data[f_new>0,])
  }
  
  return(f * sign(f)[seed_node])
}

local_density <- function(A, parameters, data)
{
  # Input: A (n x n matrix) adjacency matrix
  #        parameters (vector)
  #        -['seed node']
  #        -['connectivity']
  #        -['prune_by']
  
  # Unpack necessary parameters
  seed_node <- strtoi(parameters['seed_node'])
  connectivity <- strtoi(parameters['connectivity'])
  prune_by <- parameters['prune_by']
  n <- nrow(A)
  
  # If kNN graph, need distance to kth neighbor to prune by
  if(prune_by == 'kth_neighbor_radius')
  {
    kth_neighbor_radius <- numeric(n)
    k <- connectivity
    kNN_list <- nn2(data, k = k)$nn.idx
    kth_neighbor <- kNN_list[,k]
    for(v in 1:n)
    {
      kth_neighbor_radius[v] <- norm(data[v,], data[kth_neighbor[v],]) # distance to kth neighbor
    }
  } 
  
  # Initialize 
  distance <- rep(Inf, n)
  
  rank <- order(kth_neighbor_radius) # rank by distance to kth neighbor (smallest to largest)
  seed_rank <- which(rank %in% seed_node)
  for(v in seed_rank:n)
  {
    C <- A[rank[1:v],rank[1:v]] # pruned graph
    colnames(C) <- rank[1:v]
    components <- components(graph_from_adjacency_matrix(C))$membership # connected components of C
    seed_component <- names(which(components == components[names(components) %in% seed_node]))
    seed_component <- strtoi(seed_component) # seed connected component
    distance[seed_component] <- pmin(distance[seed_component],kth_neighbor_radius[rank[v]]) 
  }
  
  similarity <- 1 / distance
  return(similarity)
}

spectral <- function(A, parameters)
{
  # Input: A (n x n matrix) adjacency matrix,
  #        data (n x d matrix) data matrix
  # Output: second eigenvector as 'similarity' vector
  
  # unpack necessary parameters
  seed_node <- parameters['seed_node']
  
  n <- nrow(A)
  L <- laplacian(A)
  similarity <- eigen(L)$vectors[,n - 1]
  similarity <- similarity * sign(similarity[seed_node])
  return(similarity)
}

ppr <- function(A,parameters)
{
  # Input: A (n x n matrix) adjacency matrix,
  #        parameters (vector)
  #        -['seed node'] seed node
  #        -['teleportation_parameter'] teleportation parameter
  
  # unpack necessary parameters
  seed_node <- parameters['seed_node']
  alpha <- parameters['teleportation_parameter']
  
  # compute (matrix of ) ppr vector(s) (starting from each node).
  # Use formula ppr_alpha,s = 2a/(1 - a) (L + 2a/(1 - a)D)^{-1} D
  L <- laplacian(A)
  D <- diag(rowSums(A))
  similarity_matrix <- (2 * alpha / (1 - alpha)) * solve(L + 2 * alpha / (1 - alpha) * D) %*% D
  
  similarity <- similarity_matrix[seed_node,] /rowSums(A)
  
  return(similarity)
}