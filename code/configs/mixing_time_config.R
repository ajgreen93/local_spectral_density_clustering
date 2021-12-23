# sample size
n <- 20000

# parameters
sigma <- .1
gamma <- 1/10
lambda_sigma = (15/9)*sigma^gamma

# range parameters
d_range <- seq(2,3, by = 1)
D_range <- sigma*(2^seq(0,5, length.out = 10))


# family of distributions
model <- list(
  distribution_class = 'rounded_box_2_distribution',
  sigma = .1,
  gamma = gamma,
  d_range = d_range,
  D_range = D_range,
  lambda_sigma = lambda_sigma, 
  lambda = lambda_sigma * (10/9)
)

# graph configs
# connectivity should be small enough to satisfy r <= 2*sigma/sqrt(d) for all sigma and d
graph_configs <- list(
  graph_type = 'epsilon',
  connectivity = 'k_degree',
  l = model$sigma / (4 * sqrt(max(model$d_range)) )
)

empirical_configs <- list(
  n_sims = 20
)