# theoretical_pipeline.R is the main script for showing theoretical upper bounds vs. empirical samples of
# 
# Function: Iterates over sigma (thickness) and gamma (density drop) parameters
#           computing theoretical bounds and empirical estimates for
#           normalized cut.



# Dependencies
# require(MASS)
# require(ggplot2)
# require(forcats)
require(igraph)
# require(rdist)
# require(gridExtra)
# require(KernSmooth)
# require(JuliaCall)
# require(gurobi)
# require(dplyr)
# julia_setup('C:/Users/alden/AppData/Local/Julia-1.1.0/bin')
# julia_source("C:/Users/alden/OneDrive/Documents/Spectral Sparsification Work/apr.jl")

# Source methods
source_path <- 'C:/Users/alden/OneDrive/Documents/Statistics/local_spectral_density_clustering/'
code_source_path <- paste(source_path, "code/", sep = '')
script_names <- c('configs/theoretical_pipeline_config.R','theoretical_methods.R','sample.R',
                  'geometric_methods.R','graph.R')
for(script_name in script_names)
{
  source(paste(code_source_path, script_name, sep = ''))
}

# compute theoretical and sample normalized cut

# unpack necessary parameters
distribution_class <- model$distribution_class
sigma_range <- model$sigma_range
gamma_range <- model$gamma_range
d <- 2 # only considering 2 dimensions for now

theoretical_ncut <- matrix(0, ncol = length(sigma_range), nrow = length(gamma_range))
for (ii in 1:length(sigma_range))
{
  sigma <- sigma_range[ii]
  for (jj in 1:length(gamma_range))
  {
    gamma <- gamma_range[jj]
    
    # normalized cut theoretical bound 
    theoretical_ncut[jj,ii] <- ncut_upper_bound(distribution_class, sigma, gamma, d, graph_configs)
  }
}

# unpack necessary parameters
D <- model$D

empirical_ncut <- matrix(0, ncol = length(sigma_range), nrow = length(gamma_range))

# sample normalized cut
for (ii in 1:length(sigma_range))
{
  sigma <- sigma_range[ii]
  for (jj in 1:length(gamma_range))
  {
    gamma <- gamma_range[jj]
    
    # distribution
    P <- list(
      distribution_class = distribution_class,
      sigma = sigma,
      gamma = gamma,
      D = D,
      d = d
    )
    
    # sample normalized cut upper bound
    empirical_ncut[jj,ii] <- ncut_empirical_bound(P, graph_configs, empirical_configs)
    print(jj)
  }
}