# sample size
n <- 640000

# parameters
sigma <- .1
gamma <- .05
lambda_sigma = (15/9)*(sigma^gamma)

# range parameters
d_range <- c(2,3)
D_range <- sigma * 2^{seq(3,6.5,by = .5)}

# family of distributions
model <- list(
  distribution_class = 'rounded_box_2_distribution',
  sigma = sigma,
  gamma = gamma,
  d_range = d_range,
  D_range = D_range,
  lambda_sigma = lambda_sigma, 
  lambda = lambda_sigma * (10/9)
)

# graph configs
# connectivity should be small enough to satisfy r <= sigma / 4d for all sigma and d
graph_configs <- list(
  graph_type = 'epsilon',
  connectivity = 'k_degree',
  l = sigma / (4 * sqrt(max(model$d_range)) )
)

empirical_configs <- list(
  n_sims = 100
)