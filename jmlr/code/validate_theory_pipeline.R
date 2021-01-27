#---------------------------------------#
# This pipeline is used to validate our 
# theoretical bounds on
# 
# (1): normalized cut of a density cluster,
# (2): conductance of a density cluster,
# (3): local spread of a density cluster,
# (4): volume of symmetric set difference between
#      PPR estimate and density cluster,
#
# in the context of nonparametric mixtures of uniforms.
#---------------------------------------#

# How noisy should the pipeline be?
verbose <- T

# Should we include PPR?
PPR <- F
  
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
configs_name <- "2d_uniform_rectangles_epsilon"
source(file.path("configs",paste0(configs_name,".R")))

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
    psi_hat <- conductance(G_candidate_cluster,verbose)
    s_hat <- local_spread(G_candidate_cluster)
    
    # Tune hyperparameters for PPR
    seed <- which.min(dist_to_point(X,seed_location))
    alpha <- tune_alpha(G_candidate_cluster,psi_hat)
    LU <- tune_LU(G_candidate_cluster)
    
    if(PPR){
      p <- ppr(seed, alpha, G)
      deg <- degree(G)
      estimated_cluster <- best_sweep_cut(G,p/deg,LU)
      Delta_hat <- volume_ssd(empirical_candidate_cluster,estimated_cluster,G)
    } else{
      p <- rep(NA,n_samples)
      estimated_cluster <- NA
      Delta_hat <- NA
    }
    
    if(verbose & PPR & d == 2)
    {
      plot(X[,1],X[,2],xlim = c(-1,1),ylim = c(-1,1), main = c("Samples."))
      plot(X[estimated_cluster,1],X[estimated_cluster,2],xlim = c(-1,1),ylim = c(-1,1), main = c("Cluster estimate."))
      plot(X[p > 0,1],X[p > 0,2],xlim = c(-1,1),ylim = c(-1,1), main = c("Seed connected component."))
    }
    
    # Compute theoretical bounds
    phi <- ncut_bound(r,candidate_cluster,distribution)
    psi <- conductance_bound(r,candidate_cluster,distribution)
    s <- local_spread_bound(r,candidate_cluster,distribution)
    kappa <- 
      2 * phi * (8/psi^2 * log(4/s,base = exp(1)) * log(8/s,base = 2) + 4*log(8/s,base = 2) + 1)

    # Store data
    if(iter == 1)
    {
      Xs[[jj]] <- X
      ps[jj,] <- p
      empirical_candidate_clusters[jj,] <- (1:n_samples) %in% empirical_candidate_cluster
      estimated_clusters[jj,] <- (1:n_samples) %in% estimated_cluster
    }
    
    hyperparameters[[jj]][[iter]] <- list(seed = seed,
                                          alpha = alpha,
                                          LU = LU)
    empirical_ncuts[jj,iter] <- phi_hat
    empirical_conductances[jj,iter] <- psi_hat
    empirical_local_spreads[jj,iter] <- s_hat
    empirical_normalized_volume_ssd[jj,iter] <- Delta_hat / volume(G[empirical_candidate_cluster,])
    
    population_ncuts[jj] <- phi
    population_conductances[jj] <- psi
    population_local_spreads[jj] <- s
    population_condition_numbers[jj] <- kappa
    
    logger::log_info("Completed distribution ",jj,"/",n_distributions,".")
  }
  logger::log_info("Completed iteration ",iter,"/",n_iters,".")
}

# Save data
save_file = file.path("data",paste0(configs_name,".rda"))
save(list = ls(), file = save_file, compress = "xz")