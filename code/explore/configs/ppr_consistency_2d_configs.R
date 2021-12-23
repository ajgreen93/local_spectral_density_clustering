### General configs. ###
n_iters <- 5
d <- 2
ns <- round( seq(1000,4000,length.out = 10) )

### Configs for sampling. ###
n_mixtures <- 2
aspect_ratio <- 2
class_ratio <- .2
sample_X <- make_sample_mixture_of_uniforms(d, n_mixtures, aspect_ratio, class_ratio)

### Configs for neighborhood graph. ###
rs <- seq(.02,.1,by = .02)

### Methods ###
methods <- list(
  ppr = make_ppr_cluster
)

initialize_thetas <- list(ppr = make_initialize_ppr_thetas(rs))

### Error metric
error_metric <- volume_ssd
error_bounds <- list(ppr = ppr_empirical_error_bound)


