# Dependencies
library(dplyr)
library(ggplot2)
library(gridExtra)
library(magrittr)
library(Matrix)
library(RANN)
library(reshape2)
library(tibble)
library(tidyr)
source("sample.R")
source("graph.R")
source("cluster.R")
source("bound.R")
source("evaluate.R")

# Configuration
source("configs/ppr_consistency_2d_configs.R")


# Set up data structures
Xs <- vector(mode = "list", length = length(ns))
estimates <- vector(mode = "list", length = length(ns))
thetas <- vector(mode = "list", length = length(ns))
errors <- vector(mode = "list", length = length(ns))
bounds <- vector(mode = "list", length = length(ns))

for(ii in 1:length(ns))
{
  n <- ns[ii]
  
  thetas[[ii]] <- vector(mode = "list", length = length(methods))
  estimates[[ii]] <- vector(mode = "list", length = length(methods))
  bounds[[ii]] <- vector(mode = "list", length = length(methods))
  errors[[ii]] <- vector(mode = "list",length = length(methods))
  for(jj in 1:length(methods))
  {
    thetas[[ii]][[jj]] <- initialize_thetas[[jj]](sample_X,n)
    n_thetas <- nrow(thetas[[ii]][[jj]])
    estimates[[ii]][[jj]] <- vector(mode = "list",length = n_thetas)
    errors[[ii]][[jj]] <- matrix(nrow = n_thetas, ncol = n_iters)
    bounds[[ii]][[jj]] <- matrix(nrow = n_thetas, ncol = n_iters)
  }
  logger::log_info("Completed hyperparameter setup for n = ", n, ".")
}

# Estimate clusters
for(iter in 1:n_iters)
{
  save_run <- (iter == 1)
  for(ii in 1:length(ns))
  {
    n <- ns[ii]
    X <- sample_X(n)
    
    # REPLACE BY
    # empirical_cluster <- get_empirical_cluster(X,sample_X)
    cluster_center <- -c(.5,0)
    cluster_proportions <- .25*c(1,aspect_ratio)
    empirical_cluster <- get_empirical_cluster(X,cluster_center,cluster_proportions)
    
    seed <- which.min(dist_to_point(X,cluster_center))
    
    for(jj in 1:length(methods))
    {
      method <- methods[[jj]]
      bound <- error_bounds[[jj]]
      thetas_jj <- thetas[[ii]][[jj]]
      for(kk in 1:nrow(thetas_jj))
      {
        theta <- thetas_jj[kk,]
        algorithm <- method(theta, seed, empirical_cluster, error_metric)
        estimated_cluster <- algorithm(X)
        
        r <- theta[["r"]]
        G <- neighborhood_graph(X,r)
        errors[[ii]][[jj]][kk,iter] <- error_metric(estimated_cluster,empirical_cluster,G)[length(estimated_cluster)]
        bounds[[ii]][[jj]][kk,iter] <- bound(theta, seed, empirical_cluster, error_metric, G)
        if(save_run)
        {
          estimates[[ii]][[jj]][[kk]] <-  estimated_cluster
        }
      }
    }
    
    if(save_run)
    {
      Xs[[ii]] <- X
    }
    logger::log_info("Completed n = ", n, " for iter = ", iter, ".")
  }
  logger::log_info("Completed iteration ", iter, " out of ", n_iters, ".")
}

## plot errors...
path <- file.path("C:/Users/alden/OneDrive/Documents/Statistics",
  "local_spectral_density_clustering",
  "notes",
  "weekly_notes",
  "figures",
  "11_17_2020")
## ...as a function of n
error_df <- as.data.frame(matrix(NA,nrow = length(ns),ncol = length(methods)*length(rs)))
names(error_df) <- paste0("r = ",round(rs,2))
sd_df <- error_df
for(ii in 1:length(ns))
{
  for(jj in 1:length(methods))
  {
    thetas_ii_jj <- thetas[[ii]][[jj]]
    for(kk in 1:length(rs))
    {
      errors_ii_jj <- errors[[ii]][[jj]]
      errors_ii_jj_kk <- errors_ii_jj[thetas_ii_jj$r == rs[kk],,drop = F]
      min_idx <- which.min(replace_na(rowMeans(errors_ii_jj_kk),Inf))
      error_df[ii,(jj- 1)*length(rs) + kk] <- rowMeans(errors_ii_jj_kk)[min_idx]
      sd_df[ii,(jj- 1)*length(rs) + kk] <- apply(errors_ii_jj_kk,1,sd)[min_idx] / sqrt(n_iters)
    }
  }
}
plot_df <- left_join(error_df %>% 
  add_column(n = ns,.before = 1) %>%
  pivot_longer(!n, names_to = "method",values_to = "error"),
  sd_df %>% 
    add_column(n = ns,.before = 1) %>% 
    pivot_longer(!n, names_to = "method",values_to = "sd"),
  by = c("n","method"))
ggplot(plot_df,aes(x = n, y = error, colour = method)) +
  geom_line() + 
  geom_point() + 
  geom_errorbar(aes(ymin = error - sd, ymax = error + sd),width = (max(ns) - min(ns))/100) + 
  lims(y = c(0,1)) +
  labs(title = "PPR error") + 
  theme_bw() +
  theme(text = element_text(size=20))
ggsave(filename = "error.pdf",path = path,device = "pdf")

#...upper bounded
bound_df <- as.data.frame(matrix(NA,nrow = length(ns),ncol = length(methods)*length(rs)))
names(bound_df) <- paste0("r = ",round(rs,2))
for(ii in 1:length(ns))
{
  for(jj in 1:length(methods))
  {
    thetas_ii_jj <- thetas[[ii]][[jj]]
    for(kk in 1:length(rs))
    {
      bounds_ii_jj <- bounds[[ii]][[jj]]
      bounds_ii_jj_kk <- bounds_ii_jj[thetas_ii_jj$r == rs[kk],,drop = F]
      min_idx <- which.min(replace_na(rowMeans(bounds_ii_jj_kk),Inf))
      bound_df[ii,(jj- 1)*length(rs) + kk] <- rowMeans(bounds_ii_jj_kk)[min_idx]
    }
  }
}
plot_df <- bound_df %>% 
  add_column(n = ns,.before = 1) %>%
  pivot_longer(!n, names_to = "method",values_to = "bound")
ggplot(plot_df,aes(x = n, y = bound, colour = method)) +
  geom_line() + 
  geom_point() + 
  labs(title = "Allen-Zhu upper bound on PPR error") +
  theme_bw() + 
  theme(text = element_text(size=20)) +
  lims(y = c(0,50)) 
ggsave(filename = "error_upper_bound.pdf",path = path,device = "pdf")


## as a function of alpha
iis <- 10
for(ii in iis)
{
  error <- errors[[ii]][[names(methods) == "ppr"]]
  mean_error <- rowMeans(error)
  sd_error <- apply(error,1,FUN = sd) / sqrt(n_iters)
  thetas_ppr <- thetas[[ii]][[names(methods) == "ppr"]]
  alphas <- unique(thetas_ppr$alpha)
  
  print( thetas_ppr %>% add_column(error = mean_error, sd = sd_error) %>%
    ggplot(aes(x = log(alpha + .01),y = error,colour = as.factor(r))) +
    geom_line() +
    geom_point() +
    lims(y = c(0,.2)) +
    theme_bw() )
}
ggsave(filename = "error_with_alpha.pdf",path = path,device = "pdf")



## plot estimates
cluster_plots <- vector(mode = "list",length = length(ns))
preference_plots <- vector(mode = "list",length = length(ns))

iis <- 10
for(ii in iis)
{
  cluster_plots[[ii]] <- vector(mode = "list", length = length(methods))
  preference_plots[[ii]] <- vector(mode = "list", length = length(methods))
  n <- ns[[ii]]
  X <- Xs[[ii]]
  estimates_ii <- estimates[[ii]]
  thetas_ii <- thetas[[ii]]
  
  # get empirical cluster
  cluster_center <- c(-.5,0)
  cluster_proportions <- .25 * c(1,aspect_ratio)
  empirical_cluster <- get_empirical_cluster(X, cluster_center,cluster_proportions)
  
  # get seed
  seed <- which.min(dist_to_point(X,cluster_center))
  
  for(jj in 1:length(estimates_ii))
  {
    estimates_ii_jj <- estimates_ii[[jj]]
    thetas_ii_jj <- thetas_ii[[jj]]
    cluster_plots[[ii]][[jj]] <- vector(mode = "list",length = length(rs))
    preference_plots[[ii]][[jj]] <- vector(mode = "list",length = length(rs))
  
    for(kk in 1:length(rs))
    {
      errors_ii_jj <- errors[[ii]][[jj]]
      errors_ii_jj_kk <- errors_ii_jj[thetas_ii_jj$r == rs[kk],,drop = F]
      estimates_ii_jj_kk <- estimates_ii_jj[thetas_ii_jj$r == rs[kk]]
      thetas_ii_jj_kk <- thetas_ii_jj[thetas_ii_jj$r == rs[kk],]
      min_idx <- which.min(replace_na(rowMeans(errors_ii_jj_kk),Inf))
      
      estimate <- estimates_ii_jj_kk[[min_idx]]
      preference <- attr(estimate,"preference")
      
      # cluster plot
      title = paste0("r = ",round(rs[kk],3),", ",
                     "alpha = ", round(thetas_ii_jj_kk[min_idx,]$alpha,3), ".")
      plot_df <- data.frame(X)
      col <- case_when(
        1:n %in% intersect(estimate,empirical_cluster) ~ "correct",
        1:n %in% setdiff(estimate,empirical_cluster) ~ "false positive",
        1:n %in% setdiff(empirical_cluster,estimate)~  "false negative",
        attr(estimate,"preference") > 0 ~ "seed connected component",
        TRUE ~ "other"
      )
      plot_df[,"colour"] <- col
      cluster_plots[[ii]][[jj]][[kk]] <- ggplot(data = plot_df, aes(x = X1,y = X2,colour = colour)) + 
        geom_point(alpha = .1) +
        geom_point(data = data.frame(X1 = X[seed,1], X2 = X[seed,2], colour = "other"),
                   aes(x = X1, y = X2, colour = colour), pch = 3, cex = 3, show.legend = F) +
        theme_bw() +
        theme(legend.position = "none") + 
        labs(title = title) +
        scale_colour_manual(values = c("correct" = "green",
                                       "false positive" = "orange",
                                       "false negative" = "red",
                                       "seed connected component" = "blue",
                                       "other" = "black")) + coord_fixed()
      
      # preference plot
      plot_df <- data.frame(X)
      plot_df[,"colour"] <- pmin(preference,quantile(preference,.995))
      preference_plots[[ii]][[jj]][[kk]] <- ggplot(data = plot_df, aes(x = X1,y = X2,colour = colour)) + 
        geom_point(alpha = 1) + 
        theme_bw() +
        scale_colour_gradient(low = "azure2", high = "red") +
        theme(legend.position = "none") + 
        labs(title = title) + coord_fixed()
        
    }
    nCol <- floor(sqrt(length(cluster_plots[[ii]][[jj]])))
    cluster_plots_grid <- arrangeGrob(grobs = cluster_plots[[ii]][[jj]], ncol = nCol)
    preference_plots_grid <- arrangeGrob(grobs = preference_plots[[ii]][[jj]], ncol = nCol)
  }
}
ggsave(filename = "ppr_cluster.pdf",plot = cluster_plots_grid, path = path,device = "pdf")
ggsave(filename = "ppr_preference.pdf",plot = preference_plots_grid, path = path,device = "pdf")

