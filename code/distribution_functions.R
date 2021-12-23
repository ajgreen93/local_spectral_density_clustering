gaussian_mixture_sampler <- function(distribution_parameters, n)
{
  # unlist parameters
  muC <- distribution_parameters$muC
  muCpr <- distribution_parameters$muCpr
  SigmaC <- distribution_parameters$SigmaC
  SigmaCpr <- distribution_parameters$SigmaCpr
  pi <- distribution_parameters$pi
  
  # assign number of samples to each cluster
  N_class <- rmultinom(1,n,pi)
  
  z_1 <- mvrnorm(n = N_class[1], mu = muC, Sigma = SigmaC) # density of first mixture component
  z_2 <- mvrnorm(n = N_class[2], mu = muCpr, Sigma = SigmaCpr) # density of second mixture component
  z <- rbind(z_1, z_2) # density of the mixture
  
  return(z)
}

uniform_mixture_sampler <- function(sigma, epsilon, nC, nCpr)
{
  ### Specify parameters here
  n_points <- 100000
  d <- 2
  ### 1. Sample uniformly on the cube [0,5 \sigma] x [0,1]
  
  z <- matrix(ncol = 2, nrow = n_points)
  z[,1] <- runif(n_points, 0, 5 * sigma)
  z[,2] <- runif(n_points)
  
  ### 2. Assign equal probability to each point.
  ###    Take all points with x_1 \in [2 \sigma, 3 \sigma] and downweight by 1 - epsilon
  
  sampling_probabilities <- rep(1,nrow(z))
  S <- (z[,1] >= 2 * sigma) & (z[,1] <= 3 * sigma)
  sampling_probabilities[S] <- 1 - epsilon
  
  sigma_boundaries <- ( (z[,1] >= 1 * sigma) & (z[,1] <= 2 * sigma) ) | ( (z[,1] >= 3 * sigma) & (z[,1] <= 4 * sigma) )
  sampling_probabilities[sigma_boundaries] <- 1 - .00001 # we need this to be some constant amount arbitrarily lower than 1
  ### 3. sample points with replacement
  y <- z[sample.int(nrow(z), size = nC + nCpr, replace = TRUE, prob = sampling_probabilities),]
  
  return(y)
}

distribution_maker <- function(distribution_type)
{
  if(distribution_type == 'gaussian_mixture')
  {
    return(gaussian_mixture_sampler)
  } else if(distribution_type == 'uniform_mixture')
  {
    return(uniform_mixture_sample)
  } else {
    stop('You did not input a valid density.')
  }
}

gaussian_mixture_density <- function(distribution_parameters)
{
  # unlist parameters
  muC <- distribution_parameters$muC
  muCpr <- distribution_parameters$muCpr
  SigmaC <- distribution_parameters$SigmaC
  SigmaCpr <- distribution_parameters$SigmaCpr
  pi <- distribution_parameters$pi
  
  # compute density
  d <- length(muC)
  p <- function(x)
  {
    fC <- exp(-1/2 * t(x - muC) %*% solve(SigmaC) %*% (x - muC)) / sqrt( (2 * PI)^d * det(SigmaC) ) 
    fCpr <- exp(-1/2 * t(x - muCpr) %*% solve(SigmaCpr) %*% (x - muCpr)) / sqrt( (2 * PI)^d * det(SigmaCpr) )
    f <- pi[1] * fC + pi[2] * fCpr
    return( f )
  }
  return(p)
}