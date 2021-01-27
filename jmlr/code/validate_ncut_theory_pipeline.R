# How noisy should the pipeline be?
verbose <- F

# Dependencies
library(magrittr)
library(Matrix)
library(RANN)
library(reshape2)

source("bounds.R")
source("geometry.R")
source("graph.R")
source("sample.R")

# Read configs.
source("configs/2d_uniform_rectangles_epsilon.R")
n_iters <- 10

# Structures to save data.
Xs <- vector(mode = "list",length = length(distributions))
ps <- matrix(nrow = length(distributions),ncol = n_samples)
empirical_candidate_clusters <- matrix(nrow = length(distributions),ncol = n_samples)
estimated_clusters <-  matrix(nrow = length(distributions),ncol = n_samples)
hyperparameters <- vector(mode = "list",length = length(distributions))
for(jj in 1:n_distributions)
{
  hyperparameters[[jj]] <- vector(mode = "list",length = n_iters)
}

empirical_ncuts <- matrix(nrow = length(distributions),ncol = n_iters)
empirical_conductances <- matrix(nrow = length(distributions),ncol = n_iters)
empirical_local_spreads <- matrix(nrow = length(distributions),ncol = n_iters)
empirical_normalized_volume_ssd <- matrix(nrow = length(distributions),ncol = n_iters)

population_ncuts <- numeric(length(distributions))
population_conductances <- numeric(length(distributions))
population_local_spreads <- numeric(length(distributions))
population_condition_numbers <- numeric(length(distributions))

for(iter in 1:n_iters)
{
  for(jj in 1:length(distributions))
  {
    distribution <- distributions[[jj]]
    
    # Sample X_1,...,X_n
    X <- distribution$sampler(n_samples)
    
    # Form neighborhood graph
    G <- neighborhood_graph(X,r)
    
    # Find candidate cluster (and associated quantities)
    candidate_cluster <- distribution$partition(seed_location) 
    empirical_candidate_cluster <- which(apply(X,1,FUN = function(x){candidate_cluster(x)}))
    G_candidate_cluster <- G[empirical_candidate_cluster,empirical_candidate_cluster]
    
    # Compute empirical quantities
    phi_hat <- ncut(empirical_candidate_cluster,G)
    empirical_ncuts[jj,iter] <- phi_hat
    
    
    # Compute population quantities
    phi <- ncut_bound(r,candidate_cluster,distribution)
    population_ncuts[jj] <- phi
    
    logger::log_info("Completed distribution ",jj,"/",n_distributions,".")
  }
  logger::log_info("Completed iteration ",iter,"/",n_iters,".")
}