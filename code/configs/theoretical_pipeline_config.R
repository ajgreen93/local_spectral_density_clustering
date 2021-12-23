# sample size
n <- 20000

# family of distributions
model <- list(
  distribution_class = 'theory_distribution',
  sigma_range = seq(1/20, 1/10, length.out = 10),
  gamma_range = seq(1/10, 1/5, length.out = 10),
  d_range = seq(2,3, by = 1),
  D = .1
)

# graph configs
# connectivity should be small enough to satisfy r <= sigma / 4d for all sigma and d
graph_configs <- list(
  graph_type = 'epsilon',
  connectivity = min(model$sigma_range) / (4 * sqrt(max(model$d_range)) )
)

empirical_configs <- list(
  n_sims = 20
)
