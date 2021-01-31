#---------------------------------------------------------------------#
# Neighborhood graphs
#---------------------------------------------------------------------#

dist_to_point <- function(X,z)
{
  apply(X,1,FUN = function(x){sum((x - z)^2)})
}

knn_graph <- function(X,k)
{
  # unpack necessary parameters
  n <- nrow(X)
  d <- ncol(X)
  A <- Matrix(0, nrow = n, ncol = n)
  
  # compute distance matrix
  neighbors <- nn2(data = X, query = X, k = k + 1)$nn.idx[,-1]
  neighbors_list <- as.matrix(melt(neighbors)[,-2])
  A[neighbors_list] <- 1
  A <- pmax(A,t(A))
  
  return(A)
}

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

#----------------------------------------------------------------------#
# Graph functionals
#----------------------------------------------------------------------#

degree <- function(G)
{
  rowSums(G)
}

volume <- function(G)
{
  sum(G)
}

energy <- function(G,f)
{
  B_f <- balance(G,f)
  T_f <- TV(G,f)
  return(T_f / B_f)
}

TV <- function(G,f)
{
  B <- incidence(G)
  sum(abs(B %*% f))
}

balance <- function(G,f)
{
  deg <- degree(G)
  sum(deg * abs(f - weighted_median(f,wts = deg)))
}

ncut <- function(v,G)
{
  n <- nrow(G)
  min_vol <- min(volume(G[v,]), volume(G[-v,]))
  if(min_vol == 0)
  {
    Inf
  }else{
    cut(v,G) / min_vol
  }
}

cut <- function(v,G)
{
  n <- nrow(G)
  u <- setdiff(1:n,v)
  sum(G[v,u])
}

conductance <- function(G, verbose = TRUE)
{
  #-----------------------------------#
  # We compute conductance---minimum normalized cut---
  # by solving the following exact continuous relaxation of the problem:
  #
  # (1) min TV(f)/ ||f - med(f,deg)||_deg
  #
  # --- where med(f,deg) is the weighted median of f with weights proportional to degree,
  # and ||g||_deg = \sum_{i = 1}^{n} |g_i| deg_i is the degree-weighted \ell_1 norm ---
  # and taking the conductance to be the best sweep cut of f.
  # 
  # We solve (1) by suitably adapting the steepest descent algorithm of (Bresson et. al. 2013),
  # which solves a continuous relaxiation of the min ratio cut problem.
  # (https://arxiv.org/pdf/1302.2717.pdf).
  #
  # To address the problem of local minima, we use several random initializations, but
  # always take one initialization to be the Fiedler vector.
  #-----------------------------------#
  
  # Parameters
  n <- nrow(G)
  n_runs <- 0                     # How many random initializations?
  n_iters <- 500                  # How many passes over numerator and denominator of ncut?
  n_iters_inloop <- 150           # How many steps of gradient descent?
  g_0 <- fiedler_vector(G)        # Initialization at 2nd eigenvector.
  deg <- degree(G)
  theta <- .99
  tol <- 10e-10
  
  # Incidence matrix
  B <- incidence(G)
  
  # Save outputs
  ncuts <- numeric(n_runs)
  
  for(kk in 1:(n_runs + 1))
  {
    if(kk == 1)
    {
      h_0 <- g_0
    } else{
      p <- runif(n = n, min = 0, max = .33)
      h_0 <- g_0 * (2 * rbinom(n = n,size = 1,prob = p) - 1)
    }
    
    # If the 2nd eigenvalue is 0, then fiedler_vector() may not have converged.
    # In that case, the conductance is 0.
    if(is.null(g_0)){
      ncuts[kk] <- 0
      next
    }
    f_0 <- h_0 - weighted_median(h_0,w = deg); f <- f_0/norm(f_0,type = "2")
    for(ii in 1:n_iters)
    {
      # 0. Save the old
      f_old <- f
      
      # 1. identify element of subgradient
      v <- numeric(length(f))
      v[f > 0] = deg[f > 0]
      v[f < 0] = -deg[f < 0]
      v[f == 0] = (sum(deg[f < 0]) - sum(deg[f > 0])) / (sum(deg[f == 0])) * deg[f == 0]
      
      # 2. Move in the direction of subgradient
      c = 1 / (norm(v,type = "2"))
      g = f + c*v
      
      # 3. Approximately solve TV denoising problem
      #      h = argmin ||u||_{TV} + E(f)/2 ||u - g||^2
      #    until
      #      E(f) > E(h) + theta * E(f) * ||h - f||^2/||h - med(h)||
      #
      u <- g
      ubar <- u
      beta <- numeric(nrow(B))
      tau <- 1/max(degree(G))
      sigma <- tau
      lambda <- energy(G, f)
      gamma <- .7*lambda
      
      for(jj in 1:n_iters_inloop)
      {
        u_old <- u 
        beta <- pmax(-1,pmin(1,beta + sigma * B %*% ubar))
        utilde <- u - tau*(t(B) %*% beta)
        u <- (utilde + tau*lambda*g)/(1 + tau * lambda)
        sigma <- sigma * sqrt(1 + 2 * gamma * tau)
        ubar <- u + (u - u_old) / sqrt(1 + 2 * gamma * tau)
        tau <- tau / sqrt(1 + 2 * gamma * tau)
        
        if(
          energy(G, f) > energy(G, u) + theta * energy(G, f) * norm(u - f,type = "2")^2 / 
          sum(abs(deg*(u - weighted_median(u,deg))))
        ){
          break
        }
      }
      
      # 4. Center (according to degree weighted median) and normalize
      h <- as.numeric(u - weighted_median(u,deg)); f <- h/norm(h,type = "2")
      
      # 5. If you terminated without improving the energy, take the previous estimate
      if(energy(G,f) > energy(G,f_old)){
        f <- f_old
      }
      
      if(abs(energy(G,f) - energy(G,f_old))/energy(G,f_old) < tol){
        break
      }
    }
    
    # Find the optimal sweep cut
    opt_sweep_cut <- best_sweep_cut(G,f, c(-1,1))
    if(verbose)
    {
      logger::log_info("Outer loop iterations: ",ii,".")
      logger::log_info("Number of runs: ",kk,".")
      logger::log_info("Energy: ",round(energy(G,f),4),".")
      logger::log_info("Ncut: ", round(ncut(opt_sweep_cut,G),4),".")
      # plot(X[f > 0,1],X[f > 0,2],xlim = c(-1,1), ylim = c(-1,1), 
      #      main = paste0("Run number = ",kk, 
      #                    ".Outer loop iteration = ",ii, ". "))
    }
    ncuts[kk] <- ncut(opt_sweep_cut,G)
  }
  return(min(ncuts))
}

local_spread <- function(G)
{
  min(degree(G))^2/volume(G)
}

#--------------------------------------------------------------------------#
# Graph matrices
#--------------------------------------------------------------------------#

incidence <- function(G)
{
  #TODO: add capacity for weighted edges
  n <- nrow(G)
  m <- sum(G > 0)/2
  B <- Matrix(0, nrow = m, ncol = n)
  
  incidence_list <- which(G > 0, arr.ind = T);
  incidence_list <- incidence_list[incidence_list[,1] < incidence_list[,2],]
  incidence_list <- cbind(1:m,incidence_list)
  B[incidence_list[,c(1,2)]] <- 1
  B[incidence_list[,c(1,3)]] <- -1
  
  return(B)
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
  n <- nrow(A)
  1/2 * (A/ifelse(degree(A) > 0, degree(A),1) + diag(n))
}

normalized_Laplacian <- function(A)
{
  diag(nrow(A)) - replace_na(t( t(A/ sqrt(degree(A)) ) / sqrt(degree(A)) ),0)
}

#-------------------------------------------------------------#
# Graph spectral algorithms
#-------------------------------------------------------------#
fiedler_vector <- function(G,normalized = FALSE)
{
  if(normalized)
  {
    L <- normalized_Laplacian(G)
  }else{
    L <- Laplacian(G)
  }
  
  tryCatch(
    RSpectra:::eigs_sym(L, k = 2, which = "SM")$vectors[,1],
    warning = function(w){NULL}
  )
}

fiedler_value <- function(G,normalized = FALSE)
{
  if(normalized)
  {
    L <- normalized_Laplacian(G)
  }else{
    L <- Laplacian(G)
  }
  
  tryCatch(
    RSpectra:::eigs_sym(L, k = 2, which = "SM")$values[1],
    warning = function(w){0}
  )
}

ppr <- function(seed, alpha, G)
{
  # if alpha is zero, return the whole connected component of seed
  if(abs(alpha) < 1e-10){
    components <- igraph::graph.adjacency(G) %>% igraph::components()
    seed_cluster <- ifelse(components$membership == components$membership[seed],1,0)
    return(degree(G)*seed_cluster/volume(G))
  }
  
  # random walk matrix 
  W <- random_walk_Laplacian(G)
  
  # Solve the linear system (I - (1 - \alpha)W)^T p^T = alpha e_v^T
  n <- nrow(G)
  source <- ifelse(1:n == seed,1,0)
  A <- t(diag(n) - (1 - alpha)*W)
  ppr <- as.numeric(solve(A, alpha*source))
  return(ppr)
}

#-----------------------------------------------------------------#
# Clustering miscellanea
#-----------------------------------------------------------------#
best_sweep_cut <- function(G,f,bounds){
  #-------------------------------------------------------------#
  # A faster way to compute the best sweep cut.
  # Naively computing the best sweep cut --- by separately computing all sweep cuts
  # and taking the minimum --- costs O(n) matrix summations.
  #
  # By simply updating the candidate ncut as each additional vertex crosses
  # from one side of the bipartition to the other, the computation can be
  # reduced to O(n) vector summations.
  #-------------------------------------------------------------#
  n <- nrow(G)
  f_order <- order(f,decreasing = T)
  f_sort <- sort(f,decreasing = T)
  vol <- volume(G)
  
  indx_min <- max(max(which(f_sort > bounds[2])),1)
  indx_max <- min(min(which(f_sort < bounds[1])),n)
  if(indx_min == indx_max){
    return(f_order[1:indx_min])
  }
  
  cluster <- f_order[1:indx_min]
  cut <- numeric(indx_max - indx_min + 1)
  vol <- numeric(indx_max - indx_min + 1)
  cut[indx_min] <- cut(cluster,G)
  vol[indx_min] <- volume(G[cluster,])
  for(ii in (indx_min + 1):indx_max)
  {
    cluster <- c(cluster,f_order[ii])
    neighborhood_ii <- G[f_order[ii],]
    cut[ii] <- cut[ii - 1] + sum(neighborhood_ii) - 2*sum(neighborhood_ii[cluster])
    vol[ii] <- vol[ii - 1] + sum(neighborhood_ii)
  }
  ncut <- cut/pmin(vol,volume(G) - vol)
  return(f_order[1:which.min(ncut)])
}

best_sweep_cut_2 <- function(G,f,bounds)
{
  n <- nrow(G)
  f_order <- order(f,decreasing = T)
  f_sort <- sort(f,decreasing = T)
  min_ncut <- 1
  candidate_cluster <- NULL
  
  indx_min <- max(max(which(f_sort > bounds[2])),1)
  indx_max <- min(min(which(f_sort < bounds[1])),n)
  for(ii in indx_min:indx_max)
  {
    f_cluster_ii <- f_order[1:ii]
    ncut_ii <- ncut(f_cluster_ii,G)
    if(ncut_ii < min_ncut)
    {
      min_ncut <- ncut_ii
      candidate_cluster <- f_cluster_ii
    }
  }
  return(candidate_cluster)
}

volume_ssd <- function(C,S,G)
{
  volume_CminusS <- if (is.null(setdiff(C,S))) 0 else volume(G[setdiff(C,S),]) 
  volume_SminusC <- if (is.null(setdiff(S,C))) 0 else volume(G[setdiff(S,C),])
  return(volume_CminusS + volume_SminusC)
}

#----------------------------------------------------------------#
# Misc
#----------------------------------------------------------------#
weighted_median <- function(v,wts){
  median_index <- min( which( cumsum(wts[order(v)]) > sum(wts)/2 ) )
  sort(v)[median_index]
}