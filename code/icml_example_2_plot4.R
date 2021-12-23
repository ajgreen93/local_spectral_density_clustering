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
config_path <- paste(code_source_path, 'configs/icml_example3_plot1_config.R', sep = '')
distribution_path <- paste(code_source_path, 'sample.R', sep = '')
cluster_path <- paste(code_source_path, 'cluster_methods.R',sep = '')
graph_path <- paste0(code_source_path, "graph.R")
evaluation_path <- paste(code_source_path, 'evaluate.R', sep = '')

save_path <- paste0(code_source_path,'icml/')

source(distribution_path)
source(cluster_path)
source(graph_path)
source(evaluation_path)

# Read in parameters
source(config_path)

# Generate data
X <- distribution_sample(P, n)
    
# Cluster
cluster_data <- clusterer(X, cluster_methods)
    