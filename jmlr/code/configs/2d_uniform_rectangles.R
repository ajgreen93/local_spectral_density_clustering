# Distribution is 
#
# \mathbb{P} = (1 - \epsilon)/2 \cdot \sum_{k = 1}^{2} P(C_k) + epsilon \cdot P(\Omega)
#
# where \mathbb{P} is parameterized in terms of \epsilon,\sigma,\rho with
# -- \Omega = [-1,1]^2
# -- C_1 = [-1/2 - sigma/2,-1/2 + sigma/2] x [-rho/2,rho/2]
# -- C_2 = [1/2 - sigma/2,1/2 + sigma/2] x   [-rho/2,rho/2]

# for \sigma <= \rho <= 2, and \sigma < 1 - r, and r <= sigma/4d.

# Global parameters
n_samples <- 400
n_iters <- 3
d <- 2
n_mixtures <- 2
n_distributions <- 2

# PPR hyperparameters
seed_location <- c(-.5,0)
tune_LU <- function(G)
{
  return(c(1/(6 * volume(G)), 1/(5 * volume(G))))
}
tune_alpha <- function(G,psi = NULL)
{
  if(is.null(psi))
  {
    psi <- conductance(G)
  }
  s <- local_spread(G)
  alpha_theoretical <- 1/2 * 1 / (8/psi^2 * log(4/s,base = exp(1)) * log(8/s,base = 2) + 4*log(8/s,base = 2) + 1)
  return(10000*alpha_theoretical)
}

# Distributions
distributions <- vector(mode = "list",length = n_distributions)
epsilon <- .33
sigma <- .5
rhos <- seq(.5,2,length = n_distributions)
domain <- function(x)
{
  all(-1 <= x && x <= 1)
}
attributes(domain) <- list("d" = d, "measure" = 1)
for(jj in 1:n_distributions)
{
  rho <- rhos[jj]
  distributions[[jj]] <- list(
    type = "mixture_of_uniforms",
    cluster_sets =  make_two_2drectangular_cluster_sets(sigma,rhos[jj]),
    sampler = make_two_2drectangular_sampler(epsilon,sigma,rhos[jj],d),
    epsilon = epsilon,
    domain = domain,
    n_mixtures = n_mixtures
  )
  distributions[[jj]][["partition"]] <- make_partition(distributions[[jj]]$cluster_sets)
}

# Graph hyperparameters
 r <- sigma/(4 * d)
# k <- 10
