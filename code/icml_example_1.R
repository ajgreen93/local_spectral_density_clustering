# Import libraries
require(MASS)
require(ggplot2)
require(forcats)
require(igraph)
require(rdist)
require(gridExtra)
require(KernSmooth)
require(JuliaCall)
julia_setup('C:/Users/alden/AppData/Local/JuliaPro-0.6.2.1/Julia-0.6.2/bin') # takes about 20 seconds
julia_source("C:/Users/alden/OneDrive/Documents/Spectral Sparsification Work/apr.jl")

# Source methods
source_path <- 'C:/Users/alden/OneDrive/Documents/Spectral Sparsification Work/spectral-sparsification/'
code_source_path <- paste(source_path, "code/density_clustering/", sep = '')
config_path <- paste(code_source_path, 'configs/icml_example1_config.R', sep = '')
distribution_path <- paste(code_source_path, 'sample.R', sep = '')
cluster_path <- paste(code_source_path, 'cluster_methods.R',sep = '')
graph_path <- paste0(code_source_path, "graph.R")
evaluation_path <- paste(code_source_path, 'evaluate.R', sep = '')
density_cluster_path <- paste0(code_source_path, "density_cluster.R")
geometric_path <- paste0(code_source_path, "geometric_methods.R")
theoretical_path <- paste0(code_source_path, "theoretical_bounds.R")
# save_path <- paste(source_path, 'data/',gsub("[^[:alnum:]]", "", Sys.time()), sep = '')

source(distribution_path)
source(cluster_path)
source(graph_path)
source(evaluation_path)
source(density_cluster_path)
source(geometric_path)
source(theoretical_path)

# Read in parameters
source(config_path)

# Generate data
X <- distribution_sample(P, n)

# Find normalized cut
lambda_path <-  density_cluster_parameters$lambda_path
sigma <- density_cluster_parameters$sigma

Phi <- matrix(nrow = length(lambda_path),ncol = 2)
i = 1
for(lambda in lambda_path)
{
  # Compute empirical density cluster
  C <- density_cluster(P,X,density_cluster_parameters,lambda)
  
  # Find expansion of empirical density cluster
  Csig <- expansion(C,X,sigma)
  
  # Compute normalized cut of thickened set
  Phi_nr <- geometric_normalized_cut(Csig,P,X)
  
  # Compute theoretical upper bound of normalized cut of thickened set
  Phi_bf <- theoretical_upper_bound_normalized_cut(P, X, Csig, density_cluster_parameters, lambda)
  
  Phi[i,] <- c(Phi_nr, Phi_bf)
  print(Phi[i,])
  i = i + 1
}

