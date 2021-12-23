integrate <- function(f,g,d)
{
  #---------------#
  # Input: f function
  #        g grid over rectangular domain D
  #        d dimension
  # Output: a number
  # Function: estimate int_D f dx by sum_{x in g} f(x) * vol(grid_cell)
  
  # calculate necessary parameters
  l <- abs(sum(g[1,] - g[2,]))
  vol <- (l)^d
  
  return(sum(f(g)) * vol)
}

max_degree_bound <- function(P,n,r)
{
  # unpack necessary parameters
  D <- P$D
  sigma <- P$sigma
  lambda <- P$lambda
  lambda_sigma <- P$lambda_sigma
  gamma <- P$gamma
  d <- P$d
  
  # calculate additional parameters
  eta <- (lambda - lambda_sigma)/(sigma) # to compute density function q
  l <- (D + 4*sigma)/(100000^{1/d}) # to approximate integral
  
  # Upper bound expected degree
  # by E \leq r^d nu_d lambda / (int_R^d q(x) dx)
  # where the normalization comes from the fact that p(x) = q(x) / int_R^d q(x) dx
  
  ## 1. Average q over small cells in a grid, g_enl
  lim <- vector(mode = 'list', length = d)
  lim[[1]] <- seq(0, 5*sigma, by = l)
  for(ii in 2:d)
  {
    lim[[ii]] <-  seq(0, D + 4*sigma, by = l)
  }
  g_enl <- as.matrix(expand.grid(lim)) # grid over [-2sigma, 3sigma] x [0,D]^{d-1} supset C_2sigma
  
  # compute density function of each point
  q <- function(X)
  {
    # compute necessary parameters
    n_obs <- nrow(X)
    
    # compute distance of X to C
    C_vertices <- matrix(0, nrow = d, ncol = 2)
    C_vertices[1,] <- c(2*sigma,3*sigma)
    C_vertices[2:d,] <- c(rep(2*sigma,d-1), rep(D + 2*sigma,d - 1))
    dist <- dist_to_rectangle(X, C_vertices)
    
    # compute which expansion points in X are in
    C <- (dist == 0)
    C_sigma <- (dist > 0 & dist <= sigma)
    C_2sigma <- (dist > sigma & dist < 2*sigma)
    
    # compute density
    q_dens <- numeric(n_obs)
    q_dens[C] <- lambda 
    q_dens[C_sigma] <- lambda - dist[C_sigma]*eta
    q_dens[C_2sigma] <- lambda_sigma - (dist[C_2sigma] - sigma)^gamma
    
    return(q_dens)
  }
  
  c <- integrate(q,g_enl,d)
  
  ## 2. Compute upper bound on expected degree.
  nu_d <- pi^{d/2}/base::gamma(d/2 + 1) # volume of d dimensional unit ball
  E <- n*r^d*nu_d*lambda / c
  
  # upper bound max degree with probability = 1 - delta, delta = 1/20.
  M <- E + sqrt(2*E*(log(n) + log(20)))
  
  return(M)
}

min_r_bound <- function(P,n,k)
{
  #--------------------#
  # Input: P (distribution)
  #        n sample size
  #        k minimum degree
  #        integral_method flag for how we want to compute int q. 
  # Output: r, radius for which minimum degree is at least k with high probability
  #--------------------#
  
  # unpack necessary parameters
  D <- P$D
  sigma <- P$sigma
  lambda <- P$lambda
  lambda_sigma <- P$lambda_sigma
  gamma <- P$gamma
  d <- P$d
  
  # calculate additional parameters
  eta <- (lambda - lambda_sigma)/(sigma) # to compute density function q
  l <- (D + 4*sigma)/(10000000^{1/d}) # to approximate integral
  
  # Lower bound expected degree
  # by E \geq 6/25 r^d nu_d lambda_sigma / (int_R^d q(x) dx)
  # where the normalization comes from the fact that p(x) = q(x) / int_R^d q(x) dx
  
  ## 1. Average q over small cells in a grid, g_enl
  lim <- vector(mode = 'list', length = d)
  lim[[1]] <- seq(0, 5*sigma, by = l)
  for(ii in 2:d)
  {
    lim[[ii]] <-  seq(0, D + 4*sigma, by = l)
  }
  g_enl <- as.matrix(expand.grid(lim)) # grid over [-2sigma, 3sigma] x [0,D]^{d-1} supset C_2sigma
  
  # compute density function of each point
  q <- function(X)
  {
    # compute necessary parameters
    n_obs <- nrow(X)
    
    # compute distance of X to C
    C_vertices <- matrix(0, nrow = d, ncol = 2)
    C_vertices[1,] <- c(2*sigma,3*sigma)
    C_vertices[2:d,] <- c(rep(2*sigma,d-1), rep(D + 2*sigma,d - 1))
    dist <- dist_to_rectangle(X, C_vertices)
    
    # compute which expansion points in X are in
    C <- (dist == 0)
    C_sigma <- (dist > 0 & dist <= sigma)
    C_2sigma <- (dist > sigma & dist < 2*sigma)
    
    # compute density
    q_dens <- numeric(n_obs)
    q_dens[C] <- lambda 
    q_dens[C_sigma] <- lambda - dist[C_sigma]*eta
    q_dens[C_2sigma] <- lambda_sigma - (dist[C_2sigma] - sigma)^gamma
    
    return(q_dens)
  }
  
  c <- integrate(q,g_enl,d)
  
  ## 2. Compute lower bound on expected degree.
  nu_d <- pi^{d/2}/base::gamma(d/2 + 1) # volume of d dimensional unit ball
  E <- n*6/25*nu_d*lambda_sigma / c # E*r^d is a lower bound on expected degree
  
  # using Hoeffding's inequality, we obtain min(degree) \geq E - sqrt((log(n) + log(1/delta))*E) with prob 1 - \delta
  # plug in k for min degree, and solve for r
  # solve for r via quadratic formula
  A <- E
  B <- -sqrt(2*(log(n)+log(20))*E)
  C <- -k
  r <- ((-B + sqrt(B^2 - 4 * A * C))/(2*A))^{2/d}
  
  return(r)
}