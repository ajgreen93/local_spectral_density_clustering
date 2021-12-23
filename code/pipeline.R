# pipeline.R is the main *local* clustering script of local_spectral_density_clustering_project.
# Function: 1. Simulates data
#           2. Forms adjacency matrix encoding 'nearness' of data points according to a distance metric.
#           3. Iterates over clustering methods, and for each clustering method
#           -3a. Creates all hyperparameter possibilities
#           -3b. Calculates similarity scores and cluster for each hyperparameter possibility
#           -3c. Retains similarities and cluster for 'best' choice of hyperparameter
# Output: clustering_data (list)
#         - clustering_data$method (list)
#         - - clustering_data$method$similarity: similarity (vector) for given method
#         - - clustering_data$method$cluster: cluster (vector) for given method
#         - - clustering_data$method$evaluation: evaluation (numeric) for each method
  


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
script_names <- c('configs/small_config.R', 'sample.R', 'hyperparameter.R',
                  'similarity.R','cluster_methods.R','graph.R','geometric_methods.R','numeric_methods.R',
                  'evaluate.R','plot_methods.R')
for(script_name in script_names)
{
  source(paste(code_source_path, script_name, sep = ''))
}

# Generate data
data <- distribution_sample(P, n)

# Compute similarity graph
adjacency_matrix <- grapher(data, graph_configs)

# Generating local clusters
local_cluster_data <- vector(mode = 'list', length = length(local_cluster_methods))
names(local_cluster_data) <- names(local_cluster_methods)
for(cluster_method_name in names(local_cluster_methods))
{
  method <- local_cluster_methods[[cluster_method_name]]
  
  evaluation0 <- 1 # 1 should be the worst possible evaluation a cluster can receive
  
  # Form hyperparameters to optimize over
  parameter_combinations <- parameter_combiner(method, data,adjacency_matrix,graph_configs)
  
  for(ii in 1:nrow(parameter_combinations) )
  {
    parameters <- parameter_combinations[ii,]
     
    # Form cluster
    local_cluster <- local_clusterer(data, adjacency_matrix, method, parameters)
    similarity <- local_cluster$similarity
    cluster <- local_cluster$cluster
    
    # Evaluate (for tuning purposes)
    evaluation_configs <- method$evaluation_configs
    evaluation_data <- evaluator(evaluation_configs, cluster, adjacency_matrix, data, P)
    
    # Keep or discard
    evaluation <- evaluation_data$evaluation
    if(evaluation < evaluation0)
    {
      parameters0 <- parameters
      similarity0 <- similarity
      cluster0 <- cluster
      evaluation0 <- evaluation
    }
    print(paste0(c(ii,' out of ',nrow(parameter_combinations), ' completed.'), collapse = ''))
  }
  
  local_cluster_data[[cluster_method_name]] <- list(method = method,
                                 parameters = parameters0,
                                 similarity = similarity0,
                                 cluster = cluster0,
                                 evaluation = evaluation0)
}

# for a cluster tree, we estimate *one object* -- a cluster tree --
# then we optimize over cuts in that cluster tree

clustering_data <- vector(mode = 'list', length = length(cluster_tree_methods))
names(clustering_data) <- names(cluster_tree_methods)
for(cluster_method_name in names(cluster_tree_methods))
{
  method <- cluster_tree_methods[[cluster_method_name]]
  
  # estimate the cluster tree
  empirical_cluster_tree <- tree_clusterer(data, adjacency_matrix, method, parameters)
  
  evaluation0 <- 1 # 1 should be the worst possible evaluation a cluster can receive
  
  # Form hyperparameters to optimize over
  parameter_combinations <- parameter_combiner(method, data,adjacency_matrix,graph_configs)
  
  for(ii in 1:nrow(parameter_combinations) )
  {
    parameters <- parameter_combinations[ii,]
    
    # Form clustering
    clustering <- tree_to_clustering(empirical_cluster_tree, parameters, P)
    
    # Evaluate (for tuning purposes)
    evaluation_configs <- method$evaluation_configs
    evaluation_data <- evaluator(evaluation_configs, clustering, adjacency_matrix, data, P)
    
    # Keep or discard
    evaluation <- evaluation_data$evaluation
    if(evaluation < evaluation0)
    {
      parameters0 <- parameters
      similarity0 <- similarity
      clustering0 <- clustering
      evaluation0 <- evaluation
    }
    print(paste0(c(ii,' out of ',nrow(parameter_combinations), ' completed.'), collapse = ''))
  }
  
  clustering_data[[cluster_method_name]] <- list(method = method,
                                                    parameters = parameters0,
                                                    similarity = similarity0,
                                                    cluster = clustering0,
                                                    evaluation = evaluation0)
}

# Save directory
save_directory <- paste0(source_path,'data/',gsub("[^[:alnum:]]", "", Sys.time()))
data_directory <- paste0(save_directory,'/cluster_data/')
plot_directory <- paste0(save_directory,'/plots/')
for(directory in c(save_directory,data_directory,plot_directory))
{
  dir.create(directory)
}

# Produce graphics
plots <- plotter(local_cluster_data, adjacency_matrix, data, P, plot_directory)
plots <- clustering_plotter(clustering_data, adjacency_matrix, data, plot_directory)

# Save data
save_object <- list(local_cluster_data = local_cluster_data,
                    P = P)
save(save_object, file = paste0(data_directory,'local_cluster_data.R'))
