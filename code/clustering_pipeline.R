# Dependencies
require(MASS)
require(ggplot2)
require(forcats)
require(igraph)
require(rdist)
require(gridExtra)
require(KernSmooth)
require(JuliaCall)
require(gurobi)
require(dplyr)
julia_setup('C:/Users/alden/AppData/Local/Julia-1.1.0/bin')
julia_source("C:/Users/alden/OneDrive/Documents/Spectral Sparsification Work/apr.jl")

# Source methods
source_path <- 'C:/Users/alden/OneDrive/Documents/Statistics/local_spectral_density_clustering/'
code_source_path <- paste(source_path, "code/", sep = '')
script_names <- c('configs/clustering_config.R', 'sample.R', 'hyperparameter.R',
                  'similarity.R','cluster_methods.R','graph.R','geometric_methods.R',
                  'evaluate.R','plot_methods.R')
for(script_name in script_names)
{
  source(paste(code_source_path, script_name, sep = ''))
}

# Generate data
data <- distribution_sample(P, n)

# Generate similarity graph (we consider the graph as given to us a priori)
adjacency_matrix <- grapher(data, graph_configs)

# Generating clusters
cluster_data <- vector(mode = 'list', length = length(cluster_methods))
names(cluster_data) <- names(cluster_methods)

for(cluster_method_name in names(cluster_methods))
{
  method <- cluster_methods[[cluster_method_name]]
  
  evaluation0 <- 1 # 1 should be the worst possible evaluation a cluster can receive
  
  # Form hyperparameters to optimize over
  parameter_combinations <- parameter_combiner(method, data,adjacency_matrix,graph_configs)
  
  for(ii in 1:nrow(parameter_combinations) )
  {
    # Particular value of hyperparameters
    parameters <- parameter_combinations[ii,]
    
    # Make a clustering
    clustering_data <- clustering_maker(data, adjacency_matrix, method, parameters)
  }
}