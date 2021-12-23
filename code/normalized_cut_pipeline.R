# normalized_cut_pipeline.R is the main script for showing theoretical upper bounds vs. empirical samples of
# normalized cut
# Function: Iterates over sigma (thickness) and gamma (density drop) parameters
#           computing theoretical bounds and empirical estimates for
#           normalized cut.



# Dependencies
require(MASS)
require(ggplot2)
require(igraph)
require(RANN)
require(Matrix)
require(reshape2)

# Source methods
source_path <- 'C:/Users/alden/OneDrive/Documents/Statistics/local_spectral_density_clustering/'
code_source_path <- paste(source_path, "code/", sep = '')
script_names <- c('configs/diameter_normalized_cut_config.R','theoretical_methods.R','sample.R',
                  'geometric_methods.R','numeric_methods.R','graph.R')
for(script_name in script_names)
{
  source(paste(code_source_path, script_name, sep = ''))
}

# compute theoretical and sample normalized cut

# unpack necessary parameters
distribution_class <- model$distribution_class
sigma_range <- model$sigma_range
gamma <- model$gamma
lambda_sigma <- model$lambda_sigma
lambda <- model$lambda

# unpack range parameters
D_range <- model$D_range
d_range <- model$d_range

theoretical_ncut <- matrix(0, ncol = length(D_range), nrow = length(d_range))
for (ii in 1:length(D_range))
{
  D <- D_range[ii]
  for (jj in 1:length(d_range))
  {
    d <- d_range[jj]
    
    P <- list(
      distribution_class = distribution_class,
      sigma = sigma,
      gamma = gamma,
      D = D,
      d = d,
      lambda_sigma = lambda_sigma,
      lambda = lambda
    )
    # normalized cut theoretical bound 
    print(paste0('Diameter:',D, '; Dimension:', d))
    theoretical_ncut[jj,ii] <- ncut_upper_bound(P, graph_configs)
  }
}

# unpack necessary parameters
n_sims <- empirical_configs$n_sims

empirical_ncut <- array(0, dim = c(length(d_range),length(D_range),n_sims))

# sample normalized cut
for (ii in 1:length(D_range))
{
  D <- D_range[ii]
  for (jj in 1:length(d_range))
  {
    d <- d_range[jj]
    
    # distribution
    P <- list(
      distribution_class = distribution_class,
      sigma = sigma,
      gamma = gamma,
      D = D,
      d = d,
      lambda_sigma = lambda_sigma,
      lambda = lambda
    )
    
    # sample normalized cut upper bound
    print(paste0('Diameter:',D, '; Dimension:', d))
    empirical_ncut[jj,ii,] <- ncut_empirical_bound(P, graph_configs, empirical_configs)
  }
}

# save
save_directory <- paste0(source_path,'data/',gsub("[^[:alnum:]]", "", Sys.time()))
for(directory in c(save_directory))
{
  dir.create(directory)
}
save(empirical_ncut, file = paste0(save_directory, '/empirical_ncut.R'))
save(theoretical_ncut, file = paste0(save_directory, '/theoretical_ncut.R'))
configs <- list(n = n, model = model, 
                graph_configs = graph_configs,
                empirical_configs = empirical_configs)
save(configs, file = paste0(save_directory, '/configs.R'))
