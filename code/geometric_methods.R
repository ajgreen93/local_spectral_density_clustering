# Find the expansion of a geometric set
expansion <- function(C,X,sigma)
{
  # compute the minimum distance of X to C
  D_XC <- distance(X,as.matrix(X[C,]))
  distance_to_C <- apply(D_XC,1,FUN = function(d){min(d)})
  
  # return the sigma expansion set
  return(which(distance_to_C <= sigma))
}

# Compute the density cluster of the 'rounded_box' distribution over points X
rounded_box_cluster <- function(P,X)
{
  # unpack parameters
  D <- P$D
  sigma <- P$sigma
  gamma <- P$gamma
  d <- ncol(X)
  n <- nrow(X)
  
  vol_Bsigma <- (sigma^d * pi^{d/2} / base::gamma(d/2 + 1))
  c <- (1 / vol_Bsigma + d * sigma^{gamma}/(gamma + 1)) / 2^d
  l <- sigma / (8 * sqrt(d)) # grid length
  
  # set up full grid
  lim <- vector(mode = 'list', length = d)
  lim[[1]] <- seq(-2*sigma, 3*sigma, by = l )
  for(ii in 2:d)
  {
    lim[[ii]] <-  seq(-2*sigma, D + 2*sigma, by = l )
  }
  g <- as.matrix(expand.grid(lim))
  
  # set up grid over box
  lim <- vector(mode = 'list', length = d)
  lim[[1]] <- seq(0, sigma, by = l )
  for(ii in 2:d)
  {
    lim[[ii]] <-  seq(0, D, by = l )
  }
  g_box <- as.matrix(expand.grid(lim))
  
  # for each point in the grid, compute the density
  dist <- distance(g, g_box)
  grid_density_point <- matrix(0, ncol = ncol(dist), nrow = nrow(dist))
  one_sigma <- (dist <= sigma)
  two_sigma <- (dist > sigma) & (dist <= 2*sigma)
  grid_density_point[one_sigma] <- c
  grid_density_point[two_sigma] <-
    c - ((dist[two_sigma] - sigma)^{gamma}*sigma^{d - 1})/dist[two_sigma]^{d-1}
  
  grid_density <- rowSums(grid_density_point) * l^d * 1/(sigma * D^{d-1})
  
  # for each  point in X, compute the density
  dist <- distance(X, g)
  
  X_density <- numeric(n)
  for(ii in 1:n)
  {
    X_density[ii] <- grid_density[which.min(dist[ii,])]
  }
  
  # compute cluster
  lambda <- 1/(sigma * D^{d-1}) # using a cutoff of c, but this is a tuning parameter
  lambda <- 34
  C <- which(X_density >= lambda) 
  
  return(C)
}

# Compute the 'geometric' normalized cut
# C: indices of X
# X: data points
# density_cluster_parameters: self-explanatory
geometric_normalized_cut <- function(C, P, X)
{
  # unpack parameters
  d <- P$d
  lambda <- density_cluster_parameters$lambda
  sigma <- density_cluster_parameters$sigma
  r <- density_cluster_parameters$r
  
  # make neighborhood graph
  if(r == 'sigma/4d')
  {
    r <- sigma/(4*d)
  }
  D <- distance(X, X)
  kernel_function <- kernel_function_definer('uniform', r, d)
  A <- kernel_function(D)
  
  # return normalized cut
  return(normalized_cut(C,A))
}


normalized_cut <- function(S, A)
{
  # Compute the normalized cut
  # Input: S (|S| vector) subset of vertices in a graph G
  #        A (n x n matrix) adjacency matrix representation of the graph G
  if(length(S) == 0 || length(S) == nrow(A))
  {
    return(1)
  }
  cut <- sum(A[S, -S])
  vol <- min(sum(A[S,]), sum(A[-S,]))
  return(cut/vol)
}

small_side_cut <- function(S,A)
{
  # Compute E(S,S^c) / vol(S)
  # Input: S (|S| vector) subset of vertices in a graph G
  #        A (n x n matrix) adjacency matrix representation of the graph G
  # Output: numeric parameter
  if(length(S) == 0 || length(S) == nrow(A))
  {
    return(1)
  }
  cut <- sum(A[S, -S])
  vol <- sum(A[S,])
  return(cut/vol)
}

mixing_time <- function(A)
{
  # Compute the mixing time
  # Input: A (n x n) adjacency matrix representation of a graph G
  
  # unpack necessary parameters
  n <- nrow(A)
  tol <- 1/4
  e_v <- c(1, rep(0,n-1))
  
  # If contains any isolets, return infinity.
  if(components(graph_from_adjacency_matrix(A))$no > 1)
  {
    return(Inf)
  }
  
  # random walk matrix
  deg <- rowSums(A)
  W <- (1/deg * A + Diagonal(n))/2
  # W <- 1/deg * A # gambling that the Markov chain is irreducible...
  
  pi <- deg / sum(deg)
  q <- e_v %*% W
  
  r <- abs(q - pi)/pi
  m_r <- Inf # diagnostic purposes
  t = 1
  while(max(r) > 1/4)
  {
    t = t + 1
    q <- q %*% W
    r <- abs(q - pi)/pi
    
    # prevent an infinite loop
    if(m_r <= max(r))
    {
      stop('Random walk does not converge')
    } else{
      m_r <- max(r)
    }
    # print(max(r))
  }
  return(t)
}

# Compute whether or not distance between x and ellipse defined by y^2/A^2 + z^2/B^2 = 1 is <= sigma
in_expanded_ellipse <- function(x, A, B, sigma)
{
  r_x <- sqrt(x[1]^2 + x[2]^2)
  theta <- atan2(x[2],x[1])
  r_ellipse <- sqrt(1/(cos(theta)^2/sigma^2 + sin(theta)^2/D^2))
  if (r_x <= r_ellipse + sigma) return(TRUE) else return(FALSE)
}

# Euclidean norm between two vectors
norm <- function(vector_1, vector_2)
{
  return(sqrt(sum((vector_1 - vector_2)^2)))
}

intrinsic_to_extrinsic <- function(d,manifold_location)
{
  intrinsic_d <- length(manifold_location)
  extrinsic_location <- c(manifold_location,rep(0,d-intrinsic_d))
  return(extrinsic_location)
}

cluster_tree <- function(P,data, approx = T, density_data = NULL)
{
  # density related data
  if(is.null(density_data))
  {
    density_data <- density(P,data)
  }
  data_indices <- density_data$data_indices
  kNN_list <- density_data$kNN_list
  Kth_neighbor_radius <- density_data$Kth_neighbor_radius
  
  # unpack necessary parameters
  N <- nrow(kNN_list) - n
  K <- ncol(kNN_list)
  intrinsic_d <- ncol(P$centers)
  
  # Compute graph (in the form of adjacency list, ordered from vertex with smallest Kth neighbor 
  # radius to largest Kth neighbor radius, with lower indexed vertex first)
  rank <- order(Kth_neighbor_radius)
  M <- nrow(kNN_list) * ncol(kNN_list)
  adj_list <- matrix(ncol = 2, nrow = M)
  ii = 1
  for(jj in 1:(N + n))
  {
    v <- rank[jj]
    for(kk in 1:(K))
    {
      adj_list[ii,] <- c(min(v,kNN_list[v,kk]),max(v,kNN_list[v,kk]))
      ii = ii + 1
    }
  }
  density_rank <- K /  ( (N + n) * Kth_neighbor_radius[rank]^(intrinsic_d) ) ### density estimation based on manifold dimension
  
  
  max_data_rank <- max(which(rank %in% data_indices)) # no need to cut above the highest density of any sampled point
  min_data_rank <- min(which(rank %in% data_indices)) 
  if(approx == TRUE)
  {
    rank_values <- round(seq(min_data_rank, max_data_rank, length.out = 100))
    density_rank <- density_rank[rank_values]
  } else
  {
    rank_values = min_data_rank:max_data_rank
  }
  
  cluster_tree <- matrix(ncol = n, nrow = length(rank_values))
  for(ii in 1:length(rank_values) )
  {
    jj = rank_values[ii]
    high_density_nodes <- rank[1:jj]
    adj_list_C <- adj_list[adj_list[,1] %in% high_density_nodes & 
                                   adj_list[,2] %in% high_density_nodes,, drop = F]
    G <- graph_from_edgelist(adj_list_C,directed = F)
    U <- components(G)$membership[data_indices]
    
    # handle isolets
    m <- max(unique(U), na.rm = T)
    l <- sum(is.na(U))
    U[is.na(U)] <- (m + 1):(m + l)
    cluster_tree[ii,] <- U
  }
  
  return(list(empirical_cluster_tree = cluster_tree,
              density_rank = density_rank))
}

density <- function(P,data,N = 100000)
{
  # unpack necessary parameters
  intrinsic_d <- ncol(P$centers)
  n <- nrow(data)
  
  # full data to approximate the continuous cluster tree
  if(N > 0)
  {
    data_new <- distribution_sample(P,N)
    data_full <- rbind(data,data_new)
  } else{
    data_full <- data
  }
  data_indices <- 1:n
  
  # Find Kth neighbor radius for each point in the enlarged data
  Kth_neighbor_radius <- numeric(N + n)
  K <- floor(log(N + n)) # an arbitrary choice, should be order log(n) for consistent density estimation
  kNN_list <- nn2(data_full, k = K)$nn.idx
  Kth_neighbor <- kNN_list[,K]
  for(v in 1:(N + n))
  {
    Kth_neighbor_radius[v] <- norm(data_full[v,], data_full[Kth_neighbor[v],]) # distance to kth neighbor
  }
  
  return(list(data_indices = data_indices,
              kNN_list = kNN_list,
              Kth_neighbor_radius = Kth_neighbor_radius))
}

energy <- function(A,f)
{
  B_f <- balance(A,f)
  T_f <- TV(A,f)
  return(T_f / B_f / 2)
}

TV <- function(A,f)
{
  sum(abs(outer(f,f,'-')) * A)
}

balance <- function(A,f)
{
  sum(abs(f - median(f)))
}

I_Csig <- function(X, D, sigma, eps = NULL)
{
  #--------#
  # Input: matrix X, sigma, D (n x d matrix)
  # Output: vector indicator (n vector)
  # Function: I(x) =
  # {
  #    1, x in C + B(0, 2*sigma)
  #    0, otherwise
  # }
  # where C = [0,sigma] x [0,D]^{d - 1}

  # unpack parameters
  n <- nrow(X)
  d <- ncol(X)
  if(is.null(eps))
  {
    eps <- sigma
  }

  # calculate parameters
  indicator <- numeric(n)

  # determine whether points in X are in sigma expansion of C
  # compute distance of X to C
  C_vertices <- matrix(0, nrow = d, ncol = 2)
  C_vertices[1,] <- c(2*sigma,3*sigma)
  C_vertices[2:d,] <- c(rep(2*sigma,d-1), rep(D + 2*sigma,d - 1))
  dist <- dist_to_rectangle(X, C_vertices)
  return(which(dist <= (eps)))
}

dist_to_rectangle <- function(X, C_vertices)
{
  #---------------------#
  # Input: X (n x d matrix) data points
  #        C_vertices (d x 2 matrix) left and right endpoints defining an AABB
  # Output: distance of each point in X to the bounding box defined by C_vertices
  #---------------------#
  
  # unpack parameters
  n <- nrow(X)
  d <- nrow(C_vertices)
  dist_to_face <- matrix(0,nrow = n, ncol = d)
  for(jj in 1:d)
  {
    test_1 <- X[,jj] - C_vertices[jj,1]
    test_2 <- X[,jj] - C_vertices[jj,2]
    dist_to_face[,jj] <- ifelse(test_1 > 0 & test_2 < 0,0,pmin(abs(test_1), abs(test_2)))
  }
  return(sqrt(rowSums(dist_to_face^2)))
}

