neighborhood_graph <- function(X, r)
{
  # unpack necessary parameters
  d <- ncol(X)
  N <- nrow(X)
  
  k_max <- max(min(round(N*pi^{d/2}/(gamma(d/2 + 1))*r^d*2),N - 1),5)
  A <- knn_to_neighborhood_graph(X,r,k_max)
  return(A)
}

knn_to_neighborhood_graph <- function(X, r, k_max)
{
  #--------------------#
  # Input: X (n x d matrix)
  #        r connectivity parameter
  #        k_max maximum degree in graph
  # Output: A (m x n sparse adjacency matrix)
  # Function: compute adjacency matrix A
  # with A = (A_ij)_{i,j in 1:n}, A_ij = 1 iff ||X[i,] - X[j,]||_2 < r
  #--------------------#
  
  # unpack necessary parameters
  n <- nrow(X)
  d <- ncol(X)
  A <- Matrix(0, nrow = n, ncol = n)
  
  # compute distance matrix
  rneighbors <- nn2(data = X, query = X, searchtype = 'radius', k = k_max, radius = r)$nn.idx[,-1]
  rneighbors[rneighbors == 0] <- NA
  r_neighbors_list <- as.matrix(melt(rneighbors)[,-2])
  r_neighbors_list <- r_neighbors_list[rowSums(is.na(r_neighbors_list)) == 0,]
  A[r_neighbors_list] <- 1
  
  # if k_max was not high enough, at least one vertex in A will have exactly k_max - 1 neighbors.
  # Increase k_max and recalculate.
  while((max(rowSums(A)) == k_max - 1) & (k_max < n - 1))
  {
    A <- Matrix(0, nrow = n, ncol = n)
    k_max <- min(2*k_max,n - 1)
    rneighbors <- nn2(data = X, query = X, searchtype = 'radius', k = k_max, radius = r)$nn.idx[,-1]
    rneighbors[rneighbors == 0] <- NA
    r_neighbors_list <- as.matrix(melt(rneighbors)[,-2])
    r_neighbors_list <- r_neighbors_list[rowSums(is.na(r_neighbors_list)) == 0,]
    A[r_neighbors_list] <- 1
  }
  
  return(A)
}

Laplacian <- function(A)
{
  # unpack necessary parameters
  n <- nrow(A)
  
  # Laplacian matrix
  D <- Matrix(0, ncol = n, nrow = n, sparse = T)
  diag(D) <- rowSums(A)
  if(class(D - A) == "dsyMatrix")
  {
    L <- as(D - A,"matrix")
  } else{
    L <- as(D - A,"dgCMatrix")
  }
  return(L)
}

random_walk_Laplacian <- function(A)
{
  A/ifelse(degree(A) > 0, degree(A),1)
}

normalized_Laplacian <- function(A)
{
  diag(nrow(A)) - replace_na(t( t(A/ sqrt(degree(A)) ) / sqrt(degree(A)) ),0)
}

ncut <- function(v,G)
{
  cut(v,G) / min(volume(G[v,v]), volume(G[-v,-v]))
}

cut <- function(v,G)
{
  u <- setdiff(1:n,v)
  sum(G[v,u])
}

fiedler <- function(G)
{
  N <- normalized_Laplacian(G)
  tryCatch(
    RSpectra:::eigs_sym(N, k = 2, which = "SM", opts = list(retvec = FALSE))$values[1],
    warning = function(w){0}
  )
}

volume <- function(G)
{
  sum(G)
}

degree <- function(G)
{
  rowSums(G)
}

dist_to_point <- function(X,z)
{
  apply(X,1,FUN = function(x){sum((x - z)^2)})
}