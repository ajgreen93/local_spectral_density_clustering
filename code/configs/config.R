### Fundamental parameters

# Salience
sigma <- .3
epsilon <- .8

# Level set threshold
tau <- 1 / (2 * sigma * (3 - epsilon))

# Neighborhood graph radius
r <- sigma/1.01

# which cluster we will begin in
cluster_original <- 1


### Density function parameters

distribution_type <- 'gaussian_mixture' # distribution type
distribution_parameters <- list(
  muC = c(0,rep(0,1)), # mean of first mixture component
  muCpr = c(.325,rep(0,1)), # mean of second mixture component
  SigmaC = diag(2) - 
           c(1, rep(0, 1)) %*% t(c(1, rep(0, 1))) + 
           c(.08, rep(0,1)) %*% t(c(.08, rep(0, 1))), # covariance matrix of first mixture component
  SigmaCpr = diag(2) - 
    c(1, rep(0, 1)) %*% t(c(1, rep(0, 1))) + 
    c(.08, rep(0,1)) %*% t(c(.08, rep(0, 1))), # covariance matrix of second mixture component
  pi = c(.5,.5) # mixture weight distribution
)

n = 2000 # sample size
N = 6000 # number of samples needed to approximate full distribution

# Page-rank parameter grid
local_cluster_parameter_grid <- list(
  r = seq(.03,.07, by = .005), # radius for the neighborhood graph
  v = 'None',
  alpha = c(.00001, .00005, .0001, .0005, .001),
  tau = .5, # needed to find estimated density cluster, the medioid of which is used to initialize v
  cluster_size = c('.5', '.8', '1' , '1.2', '1.5', '2'),
  tol = 'None',
  normalize = TRUE
)

# Density cluster parameter grid
density_cluster_parameter_grid <- list(
  r = seq(.03,.07, by = .005),
  tau = .5,
  k = c('5','10','20','40','80')
)

true_cluster_parameters <- list(
  tau = .5
)

## simulation parameters
n_grid <- 20
n_sims <- 25