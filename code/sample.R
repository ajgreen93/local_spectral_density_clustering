# wrapper for sampling from arbitrary distribution
distribution_sample <- function(P, n)
{
  # unpack necessary parameters
  distribution_class <- P$distribution_class
  
  # use correct model
  if(distribution_class == 'gaussian')
  {
    data <- gaussian_sample(P,n)
  } else if(distribution_class == 'gaussian_mixtures')
  {
    data <- gaussian_mixture_sample(P,n)
  } else if(distribution_class == 'two_moons')
  {
    data <- two_moons_sample(P,n)
  } else if(distribution_class == 'theory_distribution')
  {
    data <- theory_sample(P,n)
  } else if(distribution_class == 'ellipse_distribution')
  {
    data <- ellipse_sample(P,n)
  } else if(distribution_class == 'rounded_box_distribution')
  {
    data <- rounded_box_sample(P,n)
  } else if(distribution_class == 'rounded_box_2_distribution')
  {
    data <- rounded_box_2_sample(P,n)
  }
  return(data)
}

rounded_box_2_sample <- function(P,n, method = 'rejection')
{
  #--------------#
  # Sample from pdf p:[0,1]^d -> (0,infty)
  # Given parameters D, sigma, lambda, lambda_sigma, gamma, d.
  # Let C = [2*sigma,3*sigma] x [2*sigma,D + 2*sigma]^{d-1}, eta = log(lambda_sigma)/log(sigma),
  #     C_sigma = C + B(0,sigma), C_2sigma = C + B(0,2*sigma).
  # q(x) = 
  # {
  #   lambda,                                x in C, 
  #   lambda - dist(x, C)^eta,               x in C_sigma \ C,
  #   lambda - sigma^eta - dist(x, C_sigma), x in C_2sigma \ C_sigma
  # }
  # Let a = (1 - int(C_2sigma) p(x)) / vol([0,1]^d / C_2sigma).
  # p(x) = 
  # {
  #   q(x), x in C_2sigma
  #   a,    x in [0,1]^d / C_2sigma
  # }
  #--------------#
  
  # unpack necessary parameters
  D <- P$D
  sigma <- P$sigma
  lambda <- P$lambda
  lambda_sigma <- P$lambda_sigma
  gamma <- P$gamma
  d <- P$d
  l <- P$l
  
  # calculate additional parameters
  eta <- (lambda - lambda_sigma)/(sigma)
  
  # compute density function of each point
  q <- function(X)
  {
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
    q_dens <- numeric(n)
    q_dens[C] <- lambda 
    q_dens[C_sigma] <- lambda - dist[C_sigma]*eta
    q_dens[C_2sigma] <- lambda_sigma - (dist[C_2sigma] - sigma)^gamma
    
    return(q_dens)
  }
  
  
  #--------------#
  # Method 1: Rejection sampling.
  #--------------#
  
  if(method == 'rejection')
  {
    
    X <- matrix(0,nrow = n, ncol = d)
    
    # sample many points within a region of interest
    # C_proposal <- [0,5*sigma] x [0, D + 4*sigma]
    
    # indexing variables
    row_counter <- 1
    sufficient_n <- FALSE
    while(sufficient_n  == FALSE)
    {
      n_proposal <- n
      X_proposal <- matrix(nrow = n_proposal, ncol = d)
      X_proposal[,1] <- runif(n_proposal,0,5*sigma)
      for(jj in 2:d)
      {
        X_proposal[,jj] <- runif(n_proposal,0,D + 4*sigma)
      }
      q_dens <- q(X_proposal)
      
      # accept / reject
      M <- lambda / (5*sigma + (D + 4*sigma)^{d - 1}) # the denominator is the uniform density over C
      P_accept <- q_dens/(M * (5*sigma + (D + 4*sigma)^{d - 1}))
      U <- runif(n_proposal,0,1)
      accept <- U <= P_accept
      X_accepted <- X_proposal[accept,]
      if(sum(accept) > n - row_counter)
      {
        X[row_counter:n,] <- X_accepted[1:(n - row_counter + 1),]
        sufficient_n <- TRUE
      } else{
        n_accept <- sum(accept)
        X[row_counter:(row_counter + n_accept - 1),] <- X_accepted
        row_counter <- row_counter + n_accept
      }
    }
    return(X)
  }
  
  # Sampling method: discretize space into cells, sample cells using multinomial, add uniform noise
  # within each cell.
  
  if(method == 'discretize')
  {
    D <- P$D
    sigma <- P$sigma
    lambda <- P$lambda
    lambda_sigma <- P$lambda_sigma
    gamma <- P$gamma
    d <- P$d
    l <- P$l
    
    # calculate additional parameters
    eta <- log(lambda - lambda_sigma)/log(sigma)
    
    # we average q over small cells in a grid, g_enl
    lim <- vector(mode = 'list', length = d)
    lim[[1]] <- seq(0, 5*sigma, by = l )
    for(ii in 2:d)
    {
      lim[[ii]] <-  seq(0, D + 4*sigma, by = l )
    }
    g_enl <- as.matrix(expand.grid(lim)) # grid over [-2sigma, 3sigma] x [0,D]^{d-1} supset C_2sigma
    
    q <- function(X)
    {
      #-------------------#
      # Input: matrix X (n x d matrix)
      # Output: vector q_dens (n vector)
      # Function: q as described above
      #-------------------#
      
      # unpack parameters
      n <- nrow(X)
      d <- ncol(X)
      
      # calculate parameters
      q_dens <- numeric(n)
      
      # set up grid over C
      lim <- vector(mode = 'list', length = d)
      lim[[1]] <- seq(2*sigma, 3*sigma, by = l )
      for(ii in 2:d)
      {
        lim[[ii]] <-  seq(2*sigma, D + 2*sigma, by = l )
      }
      g <- as.matrix(expand.grid(lim)) 
      
      # determine whether points in X are in various expansions of C
      dist <- set_distance(g, X)
      
      C <- (dist <= sqrt(d)*(l/2)) # dist <= d*(l/2)^2 implies x within C
      C_sigma <- (dist > sqrt(d)*(l/2)) & (dist <= sigma)
      C_2sigma <- (dist > sigma) & (dist <= 2*(sigma))
      
      # compute density based on distance of X to C
      q_dens[C] <- lambda 
      q_dens[C_sigma] <- lambda - dist[C_sigma]^eta
      q_dens[C_2sigma] <- lambda_sigma - (dist[C_2sigma] - sigma)^gamma
      
      return(q_dens)
    }
    
    # Now, we sample exclusively from q, rather than sampling from p. This change
    # is basically irrelevant to the theory developed, as points outside of 
    # supp(q) will not factor in to normalized cut or mixing time.
    
    # discretized cell probabilities
    cell_probs <- q(g_enl)
    cell_counts <- rmultinom(1, n, cell_probs)
    
    # cell samples
    X <- matrix(nrow = n, ncol = d)
    counter <- 1
    for(ii in 1:nrow(g_enl))
    {
      cell_count <- cell_counts[ii]
      if(cell_count > 0)
      {
        X_new <- matrix(rep(g_enl[ii,], cell_count), ncol = d, byrow = T)
        X[counter:(counter + cell_count - 1),] <- X_new
      }
      counter <- counter + cell_count
    }
    
    # add noise within each cell
    U <- runif(n*d, -l/2,l/2)
    X <- X + U
    
    return(X)
  }
}

rounded_box_sample <- function(P,n)
{
  # unpack necessary parameters
  gamma <- P$gamma
  sigma <- P$sigma
  D <- P$D
  d <- P$d
  
  X <- matrix(ncol = d, nrow = n)
  
  # 1. Sample Z
  Z <- matrix(ncol = d, nrow = n)
  Z[,1] <- runif(n,0,sigma)
  if(d >= 2)
  {
    for(ii in 2:d)
    {
      Z[,ii] <- runif(n,0,D)
    }
  }
  
  # 2. Sample epsilon
  w <- mvrnorm(n, mu = rep(0,d), Sigma = diag(d))
  
  U <- runif(n,0,1)
  r_eps <- numeric(n)
  c <- (1 + (sigma^d*sigma^gamma)/(gamma + 1)) / (sigma^d *(d + 1) / d)
  if(c < sigma^{d + gamma - 1}/(gamma + 1) )
  {
    print('Illegal density.')
  }
  p_1 <- c/d * sigma^d
  p_2 <- 1 - p_1
  L <- rbinom(n,size = 1, prob = p_1) + 1
  
  IF_1 <- function(x){c/d*x^d}
  IF_2 <- function(x){c*sigma^{d-1}*(x - sigma) - sigma^{d-1}*(x - sigma)^{gamma + 1}/(gamma + 1)}
  
  F_1 <- function(x,u){IF_1(x)/p_1 - u}
  F_2 <- function(x,u){IF_2(x)/p_2 - u}
  
  I_1 <- c(0, sigma)
  I_2 <- c(sigma, 2 * sigma)
  
  for(ii in 1:n)
  {
    if(L[ii] == 1)
    {
      r_eps[ii] <- uniroot(F_1, I_1, u = U[ii])$root
    } else
    {
      r_eps[ii] <- uniroot(F_2, I_2, u = U[ii])$root
    }
  }
  
  eps <- w / sqrt(rowSums(w^2)) * r_eps
  
  # 3. X = Z + eps
  X <- Z + eps
}

ellipse_sample <- function(P,n)
{
  # sample from a log-concave density over an expanded ellipse
  
  # unpack necessary parameters
  gamma <- P$gamma
  sigma <- P$sigma
  D <- P$D
  d <- P$d
  
  X <- matrix(ncol = d, nrow = n)
  U <- runif(n,0,1)
  c <- 1
  
  A <- (sigma/D + sigma^{gamma + 2}*(2/(gamma + 1) + 1/(gamma + 2)) + c*(sigma^3/3 + 4*sigma^2))/(9/2*sigma^2)
  for(ii in 1:nrow(X))
  {
    theta <- runif(1,0,2*pi)
  
    # region specific unscaled CDF functions
    IF_1 <- function(r){A*D*r^2/(2*sigma) - D*r^3*c/(3*sigma)}
    IF_2 <- function(r){D*r^2*(A - c)/(2*sigma)}
    IF_3 <- function(r){D*r^2*(A - c)/(2*sigma) - D/sigma * (r - 2*sigma)^{gamma + 2}/(gamma + 2) - 2*D*(r - 2*sigma)^{gamma + 1}/(gamma + 1)}
    # IF_4 <- function(r){r*(A - c - sigma^gamma)}
    
   
    # end point for fourth region
    # E_4 <- IB_sigma + sigma + (1 - 1/B_theta*(A - c/2) - 2*sigma*(A - c) + sigma^{gamma + 1}/(gamma + 1)) / 
    #   (A - c - sigma^gamma)
    
    
    # regions
    I_1 <- c(0,sigma)
    I_2 <- c(sigma, 2*sigma)
    I_3 <- c(2*sigma, 3*sigma)
    # I_4 <- c(IB_sigma + sigma, E_4)
    
    # region specific scaled CDF functions
    F_1 <- function(x,u){(IF_1(x) - IF_1(I_1[1]) )/( IF_1(I_1[2]) - IF_1(I_1[1]) )- u}
    F_2 <- function(x,u){(IF_2(x) - IF_2(I_2[1]) )/ (IF_2(I_2[2]) - IF_2(I_2[1]) ) - u}
    F_3 <- function(x,u){(IF_3(x) - IF_3(I_3[1]) )/ ( IF_3(I_3[2]) - IF_3(I_3[1]) ) - u}
    # F_4 <- function(x,u){(IF_4(x) - IF_4(I_4[1]) )/ ( IF_4(I_4[2]) - IF_4(I_4[1]) ) - u}
    
    # regional probabilities
    p_1 <- IF_1(I_1[2]) - IF_1(I_1[1])
    p_2 <- IF_2(I_2[2]) - IF_2(I_2[1])
    p_3 <- IF_3(I_3[2]) - IF_3(I_3[1])
    # p_4 <- max(IF_4(I_4[2]) - IF_4(I_4[1]),0)
    
    class <- which(rmultinom(1, size = 1, prob = c(p_1,p_2,p_3,p_4)) == 1)
    if(class == 1)
    {
      eps <- uniroot(F_1, I_1, u = U[ii])$root
    } else if (class == 2){
      eps <- uniroot(F_2, I_2, u = U[ii])$root
    } else if(class == 3){
      eps <- uniroot(F_3, I_3, u = U[ii])$root
    } else if(class == 4){
      eps <- uniroot(F_4, I_4, u = U[ii])$root 
    }
    X[ii,] <- c(eps*cos(theta), D/sigma*eps*sin(theta))
  }
  return(X)
  
}

theory_sample <- function(P,n)
{
  # unpack necessary parameters
  gamma <- P$gamma
  sigma <- P$sigma
  D <- P$D
  d <- P$d
  
  # Draw n multinomials with probabilities corresponding to the 4 distinct regions
  # of the theoretical density
  p_1 <- 2*sigma*(1 - sigma)
  p_2 <- sigma*(1 - 2*sigma) - sigma^{gamma + 1}/ (gamma + 1)
  p_3 <- sigma*(1 - 2*sigma)
  p_4 <- sigma*(1 - sigma)
  latent_class <- apply(rmultinom(n, 1, prob =  c(p_1,p_2,p_3,p_4)),2,FUN = function(x){which(x == 1)})
  
  # Draw uniforms
  U <- runif(n,0,1)
  
  # Form CDFs for each region
  F_1 <- function(x,u){(x - x^2/2)/p_1 - u}
  F_2 <- function(x,u){(x*(1 - 2*sigma) - 1 / (gamma + 1)*(x - 2*sigma)^{gamma + 1})/p_2 - (2*sigma)*(1 - 2*sigma)/p_2 - u}
  F_3 <- function(x,u){x / sigma - 3 - u}
  F_4 <- function(x,u){x / sigma - 4 - u}
  F_list <- list(F_1,F_2,F_3,F_4)
  
  # Form intervals for each region
  I_1 <- c(0, 2*sigma)
  I_2 <- c(2*sigma, 3*sigma)
  I_3 <- c(3*sigma, 4*sigma)
  I_4 <- c(4*sigma, 5*sigma)
  I_list <- list(I_1,I_2,I_3,I_4)
  
  # Draw variable from given distribution by inverting CDF for the region it is in.
  X <- matrix(ncol = d, nrow = n)
  for(ii in 1:n)
  {
    F_Z <- F_list[[latent_class[ii]]]
    I_Z <- I_list[[latent_class[ii]]]
    X[ii,] <- c(uniroot(F_Z, I_Z, u = U[[ii]])$root,c(runif(d - 1,0,D)))
  }
  return(X)
}

two_moons_sample <- function(P,n)
{
  # Input: P 'distribution' (list of parameters), including
  #        - P$d: dimension (integer greater than 0)
  #        - P$centers: 2xN matrix, each row xy coordinates for a center
  #        - P$r: radius (real number greater than 0)
  #        - P$sigma: variance matrix of noise will be simga^2 I_d
  #        - P$cluster_probabilities: probability of belonging to each moon
  # Output: data_matrix (n x d matrix)
  
  # unpack parameters
  d <- P$d
  centers <- P$centers
  r <- P$r
  sigma <- P$sigma
  cluster_probabilities <- P$cluster_probabilities
  theta_dist <- P$theta_dist
  
  # generate data lying on the two-moons manifold
  manifold_data <- matrix(0,ncol = d, nrow = n)
  
  # angles
  if(theta_dist == 'gaussian_mixture')
  {
    theta <- numeric(n)
    theta_gaussian <- matrix(c(rnorm(n,pi/4,pi/12),rnorm(n,3*pi/4, pi/12)), ncol = 2)
    theta_membership <- rbinom(n,1,.5) + 1
    for(ii in 1:n)
    {
      theta[ii] <- theta_gaussian[ii,theta_membership[ii]]
    }
  } else if(theta_dist == 'uniform'){
    theta <- runif(n,0,pi)
  }
  
  for(i in 1:n)
  {
    component <- which(rmultinom(n = 1, size = 1, prob = cluster_probabilities) == 1) # assign to a cluster
    manifold_data[i,1:2] <- centers[component,] + 
                         c(r*cos(theta[i]), 
                           (2*(component %% 2) - 1)*r*sin(theta[i])) # odd clusters are 'negative half-moons'; even clusters are 'positive half-moons'
  }
  
  # add noise around manifold
  data <- manifold_data + matrix(rnorm(n*d,0,sigma),ncol = d, nrow = n)
  
  return(data)
  
}

# sampler from gaussian distribution
gaussian_sample <- function(P, n)
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
  } else if(Sigma_type == 'None')
  {
    Sigma = P$Sigma
  }
  # print(Sigma)
  data <- mvrnorm(n = n, mu = mu, Sigma = Sigma)
  return(data)
}

# sample from mixtures of gaussian distribution
gaussian_mixture_sample <- function(P, n)
{
  # unpack necessary parameters
  d <- P$d
  cluster_probabilities <- P$cluster_probabilities
  mu_type <- P$mu_type
  Sigma_type <- P$Sigma_type
  Sigma_scale <- P$Sigma_scale
  if(!is.null(P$Sigma))
  {
    Sigma <- P$Sigma
  }else{
    Sigma <- 'None'
  }
  p <- length(cluster_probabilities)
  if(!is.null(P$mu_distance_scale))
  {
    mu_distance_scale <- P$mu_distance_scale
  }else
  {
    mu_distance_scale <- 1
  }
  
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
      shift <- (i-1)/(p * sqrt(d))*mu_distance_scale
    }
    if(Sigma_type == 'scaled')
    {
      scale <-  Sigma_scale*(1 / (4*p) )^2

    }
    P_components[[i]] <- list(d = d,
                         mu_type = 'shift',
                         shift = shift,
                         Sigma_type = Sigma_type,
                         Sigma = Sigma,
                         scale = scale)
  }
  # print(P_components[[1]])
  
  # Sample
  X <- matrix(0, ncol = d, nrow = n)
  for(i in 1:n)
  {
    component <- which(rmultinom(n = 1, size = 1, prob = cluster_probabilities) == 1)
    X[i,] <- gaussian_sample(P_components[[component]], 1)
  }
  return(X)
}