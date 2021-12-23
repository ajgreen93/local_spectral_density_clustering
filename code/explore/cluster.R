make_ppr_cluster <- function(theta,seed,empirical_cluster,error_metric)
{
  # Rank samples according to PPR on a neighborhood graph from seed,
  # output best sweep cut--- as measured by error_metric(sweep_cut,empirical_cluster)
  # as cluster estimate.
  
  r <- theta[["r"]]
  alpha <- theta[["alpha"]]
  
  ppr_cluster <- function(X){
    # graph
    G <- neighborhood_graph(X,r)
    
    # random walk matrix 
    W <- random_walk_Laplacian(G)
    
    if(abs(alpha) < 1e-10){
      components <- igraph::graph.adjacency(G) %>% igraph::components()
      cluster_estimate <- which(components$membership == components$membership[seed])
      attr(cluster_estimate,"preference") <- ifelse(1:n %in% cluster_estimate,1,0)
      return(cluster_estimate)
    }
    
    # compute PPR
    source <- ifelse(1:n == seed,1,0)
    B <- t(diag(n) - (1 - alpha)*W)
    ppr <- solve(B, alpha*source)
    ppr_normalized <- ifelse(is.na(ppr/degree(G)*volume(G)),source,ppr/degree(G)*volume(G))
    
    # estimate by thresholding
    cluster_estimate <- best_sweep_cut(ppr_normalized,empirical_cluster,error_metric, G)
    
    attr(cluster_estimate,"preference") <- as.numeric(ppr_normalized)
    return(cluster_estimate)
  }
}

make_initialize_ppr_thetas <- function(rs)
{
  initialize_ppr_thetas <- function(sample_X,n)
  {
    # Empirically choose tuning parameters for clustering.
    cluster_center <- c(-.5,0)
    cluster_proportions <- .25 * c(1,aspect_ratio)
    iters <- 10 # number of iterations over which to average
    
    # radii
    if(is.null(rs))
    {
      # set grid
      try_n <- 1000                       # if we know connectivity_radius for try_n...
      try_r_length <- 10                  # number of radii choices to try
      try_r_max <- min(cluster_proportions)/2 # maximum radius choice
      try_r_min <- try_r_max/10               # minimum radius choice
      rs_length <- 8                      # number of radii choices to end up with
      
      # find small rs which result in connected graph with high probability
      try_rs <- seq(try_r_min,try_r_max,length.out = try_r_length)
      r_mins <- numeric(iters)
      for(iter in 1:iters)
      {
        X <- sample_X(try_n)
        for(jj in 1:length(try_rs))
        {
          r <- try_rs[jj]
          G <- neighborhood_graph(X,r)
          if(!(igraph::graph.adjacency(G) %>% igraph::is.connected()))
          {
            next
          }
          r_mins[iter] <- r
          break
        }
        if(is.na(r_mins[iter]))
        {
          r_mins[iter] <- try_r_max
        }
      }
      r_min <- mean(r_mins)
      if(r_min < try_r_max)
      {
        rs <- seq(r_min,try_r_max,length.out = rs_length)
      } else{
        rs <- try_r_max
      }
      d <- environment(sample_X)$d
      
      # That was for try_n, so rescale for n.
      rs <- rs *(try_n/n * log(n)/log(try_n))^{1/d} 
    }
    
    # teleportation parameter
    iters <- 10
    thetas_list <- vector(mode = "list",length = length(rs))
    for(ii in 1:length(rs))
    {
      r <- rs[ii]
      alphas_ii <- numeric(iters)
      for(iter in 1:iters)
      {
        X <- sample_X(n)
        G <- neighborhood_graph(X,r)
        
        empirical_cluster <- get_empirical_cluster(X, cluster_center,cluster_proportions)
        G_cluster <- G[empirical_cluster,empirical_cluster]
        
        alphas_ii[iter] <- max(fiedler(G_cluster),0)
      }
      alpha_range <- c(median(alphas_ii), 10*quantile(alphas_ii,.9),quantile(alphas_ii,.1)/10)
      alphas <- c(seq(alpha_range[3],alpha_range[1],length.out = 3),
                  seq(alpha_range[1],alpha_range[2],length.out = 3))
      alphas <- unique(ifelse(alphas < 1e-10,0,alphas))
      thetas_list[[ii]] <- expand_grid(r,alpha = alphas)
    }
    thetas <- bind_rows(thetas_list)
    thetas
  }
}


get_empirical_cluster <- function(X, cluster_center,cluster_proportions)
{
  which( apply(X,1,FUN = function(x){all(abs(x - cluster_center) <= cluster_proportions)} ))
}

### Other functions. ###
best_sweep_cut <- function(p,act,error_metric,G)
{
  p_order <- order(p,decreasing = T)
  loss <- error_metric(p_order,act,G)
  p_order[1:which.min(loss)]
}

