# returns the density cluster of minimal average distance to the origin
density_cluster <- function(P, X, density_cluster_parameters, lambda)
{
  # unpack parameters
  sigma <- density_cluster_parameters$sigma
  r <- density_cluster_parameters$r
  d <- P$d
  
  # compute density level set
  Cbb_f <- density_level_set(P, X, lambda)
  
  # compute connected components of the density level set in the r-radius neighborhood graph
  if(r == 'sigma/4d')
  {
    r <- sigma/(4*d)
  }
  D <- distance(X, X)
  kernel_function <- kernel_function_definer('uniform', r, d)
  A <- kernel_function(D)
  A_level_set <- A[Cbb_f, Cbb_f]
  density_clusters <- components(graph_from_adjacency_matrix(A_level_set))
  
  origin_node <- which.min(apply(as.matrix(X[Cbb_f,]),1,FUN = function(x){sum(x^2)}))
  density_cluster <- which(density_clusters$membership == density_clusters$membership[origin_node])
  return(Cbb_f[density_cluster])
}

density_clustering <- function(P,X,density_cluster_parameters, lambda)
{
  # unpack parameters
  sigma <- density_cluster_parameters$sigma
  r <- density_cluster_parameters$r
  d <- P$d
  
  # compute density level set
  Cbb_f <- density_level_set(P, X, lambda)
  
  # compute connected components of the density level set in the r-radius neighborhood graph
  if(r == 'sigma/4d')
  {
    r <- sigma/(4*d)
  }
  D <- distance(X, X)
  kernel_function <- kernel_function_definer('uniform', r, d)
  A <- kernel_function(D)
  A_level_set <- A[Cbb_f, Cbb_f]
  density_clusters <- components(graph_from_adjacency_matrix(A_level_set))
  return(density_clusters)
}

density_level_set <- function(P, X, lambda)
{
  # compute the density function for P
  if(P$distribution_class == 'gaussian_mixtures')
  {
    f <- gaussian_mixture_density(P)
  }
  
  # evaluate the density function over X
  f_n <- apply(X, 1, FUN = function(x){f(x)})
  
  # return the indices of those points in X which have density function greater than lambda
  return(which(f_n > lambda))
}

gaussian_mixture_density <- function(P)
{
  # unpack necessary parameters
  d <- P$d
  cluster_probabilities <- P$cluster_probabilities
  mu_type <- P$mu_type
  Sigma_type <- P$Sigma_type
  Sigma_scale <- P$Sigma_scale
  p <- length(cluster_probabilities)
  
  # parameterize each mixture component
  P_components <- list()
  mu <- numeric(p)
  if(mu_type == 'evenly_spaced')
  {
    for(i in 1:p)
    {
      mu[i] <- i/p
    }
  }
  for(i in 1:p)
  {
    # defaults
    shift <- 0
    scale <- 1
    if(mu_type == 'evenly_spaced')
    {
      shift <- (i-1)/(p * sqrt(d))
    }
    if(Sigma_type == 'scaled')
    {
      scale <-  Sigma_scale * (1 / (4*p) )^2
    }
    P_components[[i]] <- list(d = d,
                              mu_type = 'shift',
                              shift = shift,
                              Sigma_type = Sigma_type,
                              scale = scale)
  }
  # print(P_components)
  
  # find density function of each mixture component
  component_densities <- vector(mode = 'list', length = p)
  for(k in 1:length(P_components))
  {
     component_densities[[k]] <- gaussian_density(P_components[[k]])
  }
  gaussian_mixture_density <- function(x)
  {
    return(sum(cluster_probabilities * sapply(component_densities, FUN = function(density){density(x)}) ) )
  }
  return(gaussian_mixture_density)
}

gaussian_mixture_densities <- function(P)
{
  # unpack necessary parameters
  d <- P$d
  cluster_probabilities <- P$cluster_probabilities
  mu_type <- P$mu_type
  Sigma_type <- P$Sigma_type
  Sigma_scale <- P$Sigma_scale
  p <- length(cluster_probabilities)
  
  # parameterize each mixture component
  P_components <- list()
  mu <- numeric(p)
  if(mu_type == 'evenly_spaced')
  {
    for(i in 1:p)
    {
      mu[i] <- i/p
    }
  }
  for(i in 1:p)
  {
    # defaults
    shift <- 0
    scale <- 1
    if(mu_type == 'evenly_spaced')
    {
      shift <- (i-1)/(p * sqrt(d))
    }
    if(Sigma_type == 'scaled')
    {
      scale <-  Sigma_scale * (1 / (4*p) )^2
    }
    P_components[[i]] <- list(d = d,
                              mu_type = 'shift',
                              shift = shift,
                              Sigma_type = Sigma_type,
                              scale = scale)
  }
  # print(P_components)
  
  # find density function of each mixture component
  component_densities <- vector(mode = 'list', length = p)
  for(k in 1:length(P_components))
  {
    component_densities[[k]] <- gaussian_density(P_components[[k]])
  }
  
  return(component_densities)
}

gaussian_density <- function(P)
{
  # unpack necessary parameters
  d <- P$d
  mu_type <- P$mu_type
  Sigma_type <- P$Sigma_type
  
  # compute specific parameters for given dimension
  if(mu_type == 'origin')
  {
    mu <- rep(0, d)
  }else if(mu_type == 'shift')
  {
    shift <- P$shift
    mu <- rep(0,d) + shift
  }
  if(Sigma_type == 'identity')
  {
    Sigma <- diag(d)
  } else if(Sigma_type == 'scaled')
  {
    scale <- P$scale
    Sigma <- diag(d) * P$scale
  }
  
  precision_matrix <- solve(Sigma)
  determinant <- det(Sigma)
  density <- function(x)
  {
    return(exp(-1/2 * t(x - mu) %*% precision_matrix %*% (x - mu)) / sqrt((2*pi)^d * determinant))
  }
  return(density)
}

