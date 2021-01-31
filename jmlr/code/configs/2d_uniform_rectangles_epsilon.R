# -----------------------------------------------------------------------------#
# Distribution is 
#
# \mathbb{P} = (1 - \epsilon)/2 \cdot \sum_{k = 1}^{2} P(C_k) + epsilon \cdot P(\Omega \ C_1 \cup C_2)
#
# where \mathbb{P} is parameterized in terms of \epsilon,\sigma,\rho with
# -- \Omega = [-1,1]^2
# -- C_1 = [-1/2 - sigma/2,-1/2 + sigma/2] x [-rho/2,rho/2]
# -- C_2 = [1/2 - sigma/2,1/2 + sigma/2] x   [-rho/2,rho/2]

# Here \sigma <= \rho < 2 - r, \sigma < 1 - r, r <= sigma/4d, and
# 0 < epsilon <= (4 - 2*rho*sigma).
# -----------------------------------------------------------------------------#

# Global parameters
n_samples <- 10000
n_iters <- 50
d <- 2
n_mixtures <- 2
n_distributions <- 10

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

## Distributions

# Sample space
domain <- function(x)
{
  all(-1 <= x && x <= 1)
}
attributes(domain) <- list("d" = d, "measure" = 4)

# Parameters
sigma <- rep(.25,n_distributions)
rho <- rep(.5,n_distributions)
density_ratio <- exp(seq(log(.1),log(.25),length.out = n_distributions))
epsilon <- density_ratio/(2*rho*sigma/(4 - 2*rho*sigma) + density_ratio)

# Make distributions
make_distribution_jj <- function(jj){
  distribution <- list(
    type = "mixture_of_uniforms",
    cluster_sets =  make_two_2drectangular_cluster_sets(sigma[jj],rho[jj]),
    sampler = make_two_2drectangular_sampler(epsilon[jj],sigma[jj],rho[jj],d),
    epsilon = epsilon[jj],
    sigma = sigma[jj],
    rho = rho[jj],
    domain = domain,
    n_mixtures = n_mixtures
  )
  distribution[["partition"]] <- make_partition(distribution$cluster_sets)
  return(distribution)
}
distributions <- lapply(1:n_distributions, make_distribution_jj)

# Graph hyperparameters
r <- min(sigma)/(4 * d)
# k <- 10

# Checks
stopifnot(all(r == r[1]))      # Don't change r.
stopifnot(all(epsilon > 0))    # Make sure it is a distribution.
stopifnot(all(rho + 2*r <= 2)) # Make sure the cut doesn't cross the domain boundary.
stopifnot(all(sigma + r <= 1)) # Make sure the cut doesn't involve 2nd cluster.
