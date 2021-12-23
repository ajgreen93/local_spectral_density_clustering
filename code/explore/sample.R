make_sample_uniform <- function(d)
{
  sample_X <- function(n){
    X <- matrix(runif(n*d,-1,1),ncol = d)
    return(X)
  }
}

make_sample_mixture_of_uniforms <- function(d, n_mixtures, aspect_ratio, class_ratio)
{
  sample_X <- function(n){
    cluster_centers <- .5 * matrix(
      c(cos(seq(0,2*pi,length.out = n_mixtures + 1))[-n_mixtures + 1], 
        sin(seq(0,2*pi,length.out = n_mixtures + 1))[-n_mixtures + 1]),
      ncol = 2
    )
    
    cluster_proportions <- .25 * matrix(
      c(rep(1,n_mixtures),
        rep(aspect_ratio,n_mixtures)),
      ncol = 2
    )
    
    Z <- make_sample_uniform(d)(n)
    probs <- c(rep(class_ratio,n_mixtures),1) %>% (function(p){p/sum(p)})
    classes <- rmultinom(n, size = 1, probs) %>% (function(M){apply(M,2,which.max)})
    X <- matrix(nrow = n, ncol = d)
    for(ii in 1:n)
    {
      if(classes[ii] <= n_mixtures)
      {
        X[ii,] <- Z[ii,]*cluster_proportions[classes[ii],] + cluster_centers[classes[ii],]
      } else{
        X[ii,] <- Z[ii,]
      }
    }
    
    X
  }
}