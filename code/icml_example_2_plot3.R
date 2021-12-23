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
config_path <- paste(code_source_path, 'configs/icml_example2_plot3_config.R', sep = '')
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
mu_distance_scale <- P_sequence$mu_distance_scale
metric_evaluation_table <- matrix(nrow = length(mu_distance_scale), ncol = n_sims)
for(i in 1:length(mu_distance_scale))
{
  for(j in 1:n_sims)
  {
    P <- P_sequence
    P$mu_distance_scale <- P_sequence$mu_distance_scale[i]
    X <- distribution_sample(P, n)
    
    # Cluster
    cluster_data <- clusterer(X, cluster_methods)
    
    # Evaluate
    metric_evaluation <- metric_evaluator(P,X,cluster_data,density_cluster_parameters)
    metric_evaluation_table[i,j] <- metric_evaluation
    print(paste0('Metric evaluation: ',metric_evaluation, ' simulation: ', j, 
                  ' mu_distance_scale: ',P$mu_distance_scale))
  }
  
}
graphical_evaluation <- graphical_evaluator_plot3(P_sequence, metric_evaluation_table)


# Save
save(plot2, file = paste0(save_path, 'plot_2.R'))
save(metric_evaluation_table, file = paste0(save_path, 'metric_evaluation_table.R'))
