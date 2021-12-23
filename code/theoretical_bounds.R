theoretical_upper_bound_normalized_cut <- function(P, X, Csig, density_cluster_parameters, lambda)
{
  ### unpack parameters
  sigma <- density_cluster_parameters$sigma
  r <- density_cluster_parameters$r
  
  ### check separation
  
  ### compute lambda_sigma
  
  # compute the density function for P
  if(P$distribution_class == 'gaussian_mixtures')
  {
    f <- gaussian_mixture_density(P)
  }
  
  # evaluate the density function over X
  f_n <- apply(as.matrix(X[Csig,]), 1, FUN = function(x){f(x)})
  lambda_sigma <- min(f_n)
  
  ### compute c_1 for low-noise condition
  
  # compute parameters and density functions for each cluster component
  mixture_components <- gaussian_mixture_density_parameterization(P)
  
  # the 'local' cluster component
  mixture_densities <- lapply(mixture_components,FUN = function(m){m$density})
  f_mixture_n <- matrix(ncol = length(mixture_densities), nrow = length(Csig))
  for(i in 1:nrow(f_mixture_n))
  {
    x <- X[Csig[i],]
    f_mixture_n[i,] <- sapply(mixture_densities, FUN = function(d){d(x)})
  }
  
  local_component_index <- which.max(colMeans(f_mixture_n))
  local_component <- mixture_components[[local_component_index]]
  
  # positive derivative term given by the 'local' cluster component
  local_Sigma <- local_component$Sigma
  local_mu <- local_component$mu
  local_density <- local_component$density
  
  distance_from_local_mu <- ginv_normal_density(local_Sigma,lambda) + sigma
  local_derivative <- distance_from_local_mu*(local_density(local_mu + distance_from_local_mu + sigma))
  
  # negative derivative term given by the 'other' components
  other_component_indices = which(1:length(mixture_components) != local_component_index)
  component_derivatives <- numeric(length(other_component_indices))
  i = 1
  for(index in other_component_indices)
  {
    component <- mixture_components[[index]]
    Sigma <- component$Sigma
    mu <- component$mu
    density <- component$density
    
    distance_from_mu <- norm(mu,local_mu) + ginv_normal_density(local_Sigma,lambda) + sigma
    component_derivatives[i] <- distance_from_mu * density(mu + abs(norm(mu,local_mu) - distance_from_local_mu - sigma))
    i = i + 1
  }
  
  # lower bound on absolute derivative
  absolute_derivative <- local_derivative - sum(component_derivatives)
  c_1 <- absolute_derivative / 2
  
  ### return the ub of Theorem 1
  d <- P$d
  if(r == 'sigma/4d')
  {
    r <- sigma/(4*d)
  }
  return(r * 4 / sigma * d * lambda/lambda_sigma * (lambda - c_1*r/2) / lambda_sigma)
}

# compute the distance of the lambda level set to a the centroid of a normal distribution
# with standard deviation tau
# TODO: modify to work for greater than one dimension
ginv_normal_density <- function(Sigma, lambda)
{
  sqrt( 2*log(1 / (lambda * sqrt(2 * pi * Sigma) ) ) * Sigma )
}

gaussian_mixture_density_parameterization <- function(P)
{
  # unpack necessary parameters
  d <- P$d
  cluster_probabilities <- P$cluster_probabilities
  mu_type <- P$mu_type
  Sigma_type <- P$Sigma_type
  Sigma_scale <- P$Sigma_scale
  p <- length(cluster_probabilities)
  
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
      shift <- (i-1)/p
    }
    if(Sigma_type == 'scaled')
    {
      scale <-  Sigma_scale * (d / (4*p) )^2
    }
    P_components[[i]] <- list(d = d,
                              mu_type = 'shift',
                              shift = shift,
                              Sigma_type = Sigma_type,
                              scale = scale)
  }
  
  mixture_components <- vector(mode = 'list', length = p)
  for(k in 1:length(P_components))
  {
    mixture_components[[k]] <- gaussian_density_parameterization(P_components[[k]])
  }
  return(mixture_components)
}

gaussian_density_parameterization <- function(P)
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
  return(list(density = density,
              mu = mu,
              Sigma = Sigma))
}