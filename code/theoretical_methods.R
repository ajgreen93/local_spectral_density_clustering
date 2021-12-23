ncut_upper_bound <- function(P, graph_configs)
{
  # Input: - distribution_class (flag),
  #        - sigma (number),
  #        - gamma (number),
  #        - graph_configs (list)
  
  ## unpack necessary parameters
  sigma <- P$sigma
  lambda <- P$lambda
  lambda_sigma <- P$lambda_sigma
  gamma <- P$gamma
  d <- P$d
  graph_type <- graph_configs$graph_type
  r <- graph_configs$connectivity
  
  # calculate necessary parameters
  if(r == 'theoretical_max')
  {
    r <- sigma / 4*sqrt(d)
  } else if (r == 'k_degree') {
    r <- min_r_bound(P,n,1)
  }
  # enforce the condition on r required by normalized cut upper bound
  if(r >= (3*sigma) / (2*d))
  {
    r <- (3*sigma) / (2*d)
  }
  
  bound <- r * 4 / sigma * d * lambda / lambda_sigma * (lambda_sigma - r^{gamma}/(gamma + 1)) /
    lambda_sigma
  return(bound)
}

mixing_time_upper_bound <- function(P, graph_configs)
{
  # Input: - distribution_class (flag),
  #        - sigma (number),
  #        - gamma (number),
  #        - graph_configs (list)
  
  ## unpack necessary parameters
  lambda <- P$lambda
  lambda_sigma <- P$lambda_sigma
  Lambda_sigma <- lambda # TODO: change this; it holds only when the density is uniform over C.
  D <- P$D
  d <- P$d
  r <- graph_configs$connectivity
  
  # calculate necessary parameters
  if(r == 'theoretical_max')
  {
    r <- sigma / 4*sqrt(d)
  } else if (r == 'k_degree') {
    r <- min_r_bound(P,n,1)
  }
  
  # r must be less than the maximum allowed by theory
  if(r >= 4 * sigma / sqrt(d) )
  {
    r <- 4*sigma / sqrt(d)
  }
  
  lambda_ratio <- Lambda_sigma / lambda_sigma
  c_lambda <- 2*log(lambda_ratio)
  mu <- log(2*(D + 2*sigma) / r)
  
  bound <- (lambda_ratio)^8 * (mu + c_lambda)*(c_lambda*(d^3*mu + 
              d^2*(D + 2*sigma)^2/r^2*(log(d) + lambda_ratio^2 * 2 * log(lambda_ratio) + log(mu)))) # gives the bound up to multiplicative constants
  
  return(bound)
}

ncut_empirical_bound <- function(P,  graph_configs, empirical_configs)
{
  # Input: - distribution_class (flag),
  #        - sigma (number),
  #        - gamma (number),
  #        - graph_configs (list)
  
  ## unpack necessary parameters
  graph_type <- graph_configs$graph_type
  connectivity <- graph_configs$connectivity
  n_sims <- empirical_configs$n_sims
  
  ncut <- numeric(n_sims)
  for(ii in 1:n_sims)
  {
    if(distribution_class == 'theory_distribution')
    {
      # sample from the theoretical distribution
      X <- distribution_sample(P,n)
      
      # compute graph from sample
      sigma <- P$sigma
      d <- ncol(X)
      
      adjacency_matrix <- grapher(X, graph_configs, sigma)
      
      # compute normalized cut from graph
      Csig <- which(X[,1] <= 2*sigma)
      ncut[ii] <- normalized_cut(Csig,adjacency_matrix)
      print(ncut[ii])
    } else if(distribution_class == 'rounded_box_2_distribution')
    {
      # sample from the theoretical distribution
      l <- graph_configs$l
      P$l <- l
      X <- distribution_sample(P,n)
      
      # neighborhood graph
      r <- graph_configs$connectivity
      if(r == 'theoretical_max')
      {
        d <- P$d
        r <- sigma / (4 * sqrt(d))
      } else if(r == 'k_degree')
      {
        r <- min_r_bound(P,n,1)
      }
      # enforce the condition on r required by normalized cut upper bound
      if(r >= (3*sigma) / (2*d))
      {
        r <- (3*sigma) / (2*d)
      }
      M <- max(ceiling(max_degree_bound(P,n,r)),3) # needed to initialize kNN graph for large enough k
      adjacency_matrix <- sparse_neighborhood_grapher(X, r = r,k_max = M)
      
      # compute Csig
      D <- P$D
      sigma <- P$sigma
      Csig <- I_Csig(X,D,sigma)
      
      # compute normalized cut of subgraph induced by Csig
      ncut[ii] <- small_side_cut(Csig, adjacency_matrix)
      print(ncut[ii])
    }
  }
  return(ncut)
}

mixing_time_empirical_bound <- function(P, n, graph_configs, empirical_configs)
{
  # Input: - distribution_class (flag),
  #        - sigma (number),
  #        - gamma (number),
  #        - graph_configs (list)
  # Output: average number of steps needed for random walk to mix
  
  ## unpack necessary parameters
  graph_type <- graph_configs$graph_type
  connectivity <- graph_configs$connectivity
  n_sims <- empirical_configs$n_sims
  
  tau <- numeric(n_sims)
  for(ii in 1:n_sims)
  {
    if(distribution_class == 'theory_distribution')
    {
      # sample from the theoretical distribution
      X <- distribution_sample(P,n)
      
      # compute graph over Csig from sample
      sigma <- P$sigma
      d <- ncol(X)
      Csig <- which(X[,1] <= 2*sigma)
      Csig_adjacency_matrix <- grapher(X[Csig,], graph_configs, sigma)
      
      # compute mixing time from graph
      tau[ii] <- mixing_time(Csig_adjacency_matrix)
      print(tau[ii])
    }
    if(distribution_class == 'ellipse_distribution')
    {
      # sample from the theoretical distribution
      X <- distribution_sample(P,n)
      
      # compute Csig
      sigma <- P$sigma
      D <- P$D
      Csig <- which(apply(X,1,FUN = function(x){in_expanded_ellipse(x,sigma,D,sigma)}))
      
      # induced neighborhood subgraph
      Csig_adjacency_matrix <- grapher(X[Csig,], graph_configs, sigma)
      
      tau[ii] <- mixing_time(Csig_adjacency_matrix)
    } else if(distribution_class == 'rounded_box_distribution')
    {
      # sample from the theoretical distribution
      X <- distribution_sample(P,n)
      
      # compute C
      C <- rounded_box_cluster(P,X)
      
      # compute Csig
      sigma <- P$sigma
      C_sig <- expansion(C,X,sigma)
    } else if(distribution_class == 'rounded_box_2_distribution')
    {
      X <- distribution_sample(P,n)
      
      # compute Csig
      D <- P$D
      sigma <- P$sigma
      Csig <- I_Csig(X,D,sigma)
      
      # induced neighborhood subgraph
      r <- graph_configs$connectivity
      if(r == 'theoretical_max')
      {
        d <- P$d
        r <- sigma / (4 * sqrt(d))
      } else if(r == 'k_degree')
      {
        r <- min_r_bound(P,n,1)
      }
      # r must be less than the maximum allowed by theory
      if(r >= 4 * sigma / sqrt(d) )
      {
        r <- 4*sigma / sqrt(d)
      }
      M <- ceiling(max_degree_bound(P,n,r))
      Csig_adjacency_matrix <- sparse_neighborhood_grapher(X[Csig,], r = r, k_max = M)
      
      if(!isSymmetric(Csig_adjacency_matrix))
      {
        stop('The r-neighborhood graph was made incorrectly.')
      }
      
      print(Sys.time())
      tau[ii] <- mixing_time(Csig_adjacency_matrix)
      print(tau[ii])
    }
    print(paste0(ii,' out of ',n_sims, ' simulations completed.'))
  }
  return(tau)
}