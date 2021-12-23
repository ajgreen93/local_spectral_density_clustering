library(RANN)
grapher <- function(data, graph_configs, sigma = NULL)
{
  # Input: data (matrix)
  #        graph_configs (list of parameters)
  #        - graph_configs$graph_type: flag for type of graph ('kNN', 'epsilon', 'Gaussian')
  #        - graph_configs$connectivity: minimum number of neighbors each point should have
  
  # unpack parameters
  graph_type <- graph_configs$graph_type
  connectivity <- graph_configs$connectivity
  d <- ncol(data)
  
  if(graph_type == 'kNN')
  {
    adjacency_matrix <- knn_grapher(data, connectivity)
  } else if(graph_type == 'epsilon')
  {
    if(connectivity == 'theoretical_max')
    {
      # derive value of r from theory
      r <- sigma / (4 * d)
    } else{
      r <- connectivity
    }
    cluster_method <- list(
      kernel = 'uniform',
      h = r
    )
    adjacency_matrix <- kernel_grapher(data, cluster_method)
  } else if(graph_type == 'k_to_epsilon')
  {
    adjacency_matrix <- k_to_r_graph(connectivity,data)
  }
  return(adjacency_matrix)
}

knn_grapher <- function(data,k)
{
  # Input: data (n x d matrix)
  #        k: number of neighbors
  # Output: n x n adjacency matrix, with (u ~ v) if either u is a jth-nearest neighbor for v
  #                                 or v is a jth-nearest neighbor for u
  
  # compute necessary parameters
  n <- nrow(data)
  
  # form kNN graph
  kNN_list <- nn2(data, k = k)$nn.idx
  A <- matrix(0, nrow = n, ncol = n)
  for(i in 1:n)
  {
    u <- i
    N_i <- kNN_list[i,]
    for(j in 1:k)
    {
      v <- N_i[j]
      A[u,v] <- 1
      A[v,u] <- 1 # (u,v) in the graph if u is a nearest neighbor of v or vice versa
    }
  }
  return(A)
}

set_distance <- function(X,Y = NULL)
{
  #--------------------#
  # Input: X (n x d matrix)
  #        Y (m x d matrix)
  #        max_dist parameter
  # Output: dist (length m vector)
  # Function: compute distance vector dist with
  # dist_i = min_{j in 1:n} ||Y[i,] - X[j,]||_2
  
  # if no query matrix, query against itself
  if(is.null(Y))
  {
    Y <- X
  }
  
  # unpack necessary parameters
  n <- nrow(X)
  m <- nrow(Y)
  d <- ncol(X)
  A <- Matrix(0, nrow = m, ncol = n)
  
  dist <- nn2(data = X, query = Y, k = 1)$nn.dists
  
  return(dist)
}

sparse_neighborhood_grapher <- function(X, Y = NULL, r, k_max)
{
  #--------------------#
  # Input: X (n x d matrix)
  #        Y (m x d matrix)
  #        r connectivity parameter
  #        k_max maximum degree in graph
  # Output: A (m x n sparse adjacency matrix)
  # Function: compute adjacency matrix A
  # with A = (A_ij)_{i,j in 1:n}, A_ij = 1 iff ||Y[i,] - X[j,]||_2 < r
  
  # if no query matrix, query against itself
  if(is.null(Y))
  {
    Y <- X
  }
  
  # unpack necessary parameters
  n <- nrow(X)
  m <- nrow(Y)
  d <- ncol(X)
  A <- Matrix(0, nrow = m, ncol = n)
  
  # compute distance matrix
  rneighbors <- nn2(data = X, query = Y, searchtype = 'radius', k = k_max, radius = r)$nn.idx[,-1]
  rneighbors[rneighbors == 0] <- NA
  r_neighbors_list <- as.matrix(melt(rneighbors, na.rm = T)[,-2])
  A[r_neighbors_list] <- 1
  
  # if k_max was not high enough, at least one vertex in A will have exactly k_max - 1 neighbors.
  # Increase k_max and recalculate.
  while(max(rowSums(A)) == k_max - 1)
  {
    A <- Matrix(0, nrow = m, ncol = n)
    k_max <- k_max + 50
    rneighbors <- nn2(data = X, query = Y, searchtype = 'radius', k = k_max, radius = r)$nn.idx[,-1]
    rneighbors[rneighbors == 0] <- NA
    r_neighbors_list <- as.matrix(melt(rneighbors, na.rm = T)[,-2])
    A[r_neighbors_list] <- 1
  }
  
  return(A)
}

kernel_grapher <- function(X, cluster_method)
{
  # compute necessary parameters
  n <- nrow(X)
  d <- ncol(X)
  
  # form kernel graph
  D <- distance(X, X)
  if(cluster_method$h == 'quantile')
  {
    quantile <- cluster_method$quantile
    cluster_method$h <- quantile(D, quantile)
    print(cluster_method$h)
  }
  kernel_function <- kernel_function_definer(cluster_method$kernel, cluster_method$h, d)
  A <- kernel_function(D)
  diag(A) <- 0
  
  return(A)
}

distance <- function(X,Y)
{
  return( sqrt(matrix(rep(diag(Y %*% t(Y)), nrow(X)), ncol = nrow(Y), nrow = nrow(X), byrow = T) +
                 matrix(rep(diag(X %*% t(X)), nrow(Y)), ncol = nrow(Y), nrow = nrow(X)) -
                 2 *(X %*% t(Y))) )
}

# NOTE: The resulting kernel function should take in a matrix of distances
#       and return pairwise evaluations of a kernel function of the form
#       k(x,y) = k(||x - y||).
kernel_function_definer <- function(kernel, h, d)
{
  if(kernel == 'gaussian')
  {
    kernel_function <- function(D)
    {
      exp(-1/2 * D^2 / h^2)
    }
  } else if(kernel == 'uniform')
  {
    kernel_function <- function(D)
    {
      (D <= h)
    }
  }
  return(kernel_function)
}

laplacian <- function(A)
{
  D <- diag(rowSums(A))
  L <- D - A
  return(L)
}

laplacian_solver <- function(L)
{
  dnx <- dimnames(L)
  if(is.null(dnx)) dnx <- vector("list", 2)
  s <- svd(L)
  tol = sqrt(.Machine$double.eps)
  nz <- s$d > tol * s$d[1]
  nu <- 1 / s$d[nz]
  Linv <- structure(
    if(any(nz)) s$v[, nz] %*% (t(s$u[, nz]) * nu) else L,
    dimnames = dnx[2:1])
  return(list(Linv = Linv,
              nu = nu) )
}

incidence <- function(A)
{
  #TODO: add capacity for weighted edges
  B <- matrix(0, nrow = sum(A > 0), ncol = ncol(A) )
  A_list <- which(A > 0, arr.ind = T)
  A_list <- A_list[which(A_list[,1] <= A_list[,2]),]
  for(i in 1:nrow(A_list) )
  {
    edge <- A_list[i,]
    B[i,edge[1]] <- 1
    B[i,edge[2]] <- -1
  }
  return(B)
}

k_to_r_graph <- function(k,data)
{
  #--------------------#
  # Input: k, minimum number of neighbors each vertex should have
  #        data
  # Output: r, radius needed to achieve at k neighbors for each vertex
  #--------------------#
  
  # distance to kth nearest neighbor
  dist_to_k <- nn2(data, k = k)$nn.dists[,k]
  
  # r should be at least the maximum of such distances, to guarantee no isolets
  r <- sort(dist_to_k, decreasing = T)[1]
  
  # estimated maximum number of neighbors
  k_max <- n/10
  
  # adjacency matrix
  A <- sparse_neighborhood_grapher(data,Y = NULL,r,k_max)
  return(A)
}

k_to_max_k <- function(k,data)
{
  
}