# sample size
n <- 20000

# parameters
sigma <- .1
D <- sigma * 2^5
gamma <- 1/10
lambda_sigma = (15/9)*sigma^gamma

# family of distributions
model <- list(
  distribution_class = 'rounded_box_2_distribution',
  sigma_range = sigma * 2^{seq(0,5,length.out = 10)},
  gamma = gamma,
  d_range = seq(2,3, by = 1),
  D = D,
  lambda_sigma = lambda_sigma, 
  lambda = lambda_sigma * (10/9)
)

# graph configs
# connectivity should be small enough to satisfy r <= 2*sigma/sqrt(d) for all sigma and d
graph_configs <- list(
  graph_type = 'epsilon',
  connectivity = 'k_degree',
  l = min(model$sigma_range) / (4 * sqrt(max(model$d_range)) )
)

empirical_configs <- list(
  n_sims = 50
)