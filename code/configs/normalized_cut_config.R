# sample size
n <- 20000

# parameters
sigma <- .1
D <- sigma * 2^5
gamma <- .05
lambda_sigma = (15/9)*(D^gamma)
d_range <- c(2,3)


# family of distributions
model <- list(
  distribution_class = 'rounded_box_2_distribution',
  sigma_range = sigma * 2^{seq(0,5,length.out = 10)},
  gamma = .05,
  d_range = d_range,
  D = D,
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