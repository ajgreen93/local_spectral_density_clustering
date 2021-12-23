# mixing_time_pipeline.R is the main script for showing theoretical upper bounds vs. empirical
# bounds for mixing time over a density cluster
# Function: Iterates over sigma (thickness) and gamma (density drop) parameters
#           computing theoretical bounds and empirical estimates for
#           normalized cut and mixing time.



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
script_names <- c('configs/diameter_mixing_time_config.R','theoretical_methods.R','sample.R',
                  'geometric_methods.R','numeric_methods.R','graph.R')
for(script_name in script_names)
{
  source(paste(code_source_path, script_name, sep = ''))
}

# unpack necessary parameters
distribution_class <- model$distribution_class
sigma <- model$sigma
gamma <- model$gamma
lambda_sigma <- model$lambda_sigma
lambda <- model$lambda

# unpack range parameters
d_range <- model$d_range
D_range <- model$D_range

# compute theoretical mixing time 
theoretical_mixing_time <- matrix(0, ncol = length(D_range), nrow = length(d_range))
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
    
    print(paste0('Dimension:',d, '; Diameter:', D))
    # mixing time theoretical bound 
    theoretical_mixing_time[jj,ii] <- mixing_time_upper_bound(P, graph_configs)
  }
}

# compute empirical mixing time
n_sims <- empirical_configs$n_sims
empirical_mixing_time <- array(0, dim = c(length(D_range),length(d_range), n_sims) )
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
  
    print(paste0('Dimension:',d, '; Diameter:', D))
    # mixing theoretical bound 
    empirical_mixing_time[ii,jj,] <- mixing_time_empirical_bound(P, n, graph_configs, empirical_configs)
  }
}

# save
save_directory <- paste0(source_path,'data/',gsub("[^[:alnum:]]", "", Sys.time()))
for(directory in c(save_directory))
{
  dir.create(directory)
}
save(empirical_mixing_time, file = paste0(save_directory, '/empirical_mixing_time.R'))
save(theoretical_mixing_time, file = paste0(save_directory, '/theoretical_mixing_time.R'))
configs <- list(n = n, model = model, 
                graph_configs = graph_configs,
                empirical_configs = empirical_configs)
save(configs, file = paste0(save_directory, '/configs.R'))