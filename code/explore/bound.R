ppr_empirical_error_bound <- function(theta, seed, empirical_cluster, error_metric, G)
{
  alpha <- theta$alpha
  (15 * ncut(empirical_cluster, G) / alpha  + 8  * ncut(empirical_cluster,G))
}

ppr_theoretical_error_bound <- function(theta)
{
  r <- theta$r
  h <- .25*aspect_ratio*2
  w <- .25*2
  epsilon <- (c(rep(class_ratio,n_mixtures),1) %>% (function(p){p/sum(p)}))[n_mixtures+1]
  
  p_out <- epsilon / 4
  p_in <-  (epsilon / 4 + (1 - epsilon)/(h*w))
  ncut_predicted <- r * p_out / p_in * 2*(w + h)/(w*h)
  
  fiedler_predicted <- 4 * r^2/h^2 / 2
  
  condition_number <- ncut_predicted / fiedler_predicted
  
}