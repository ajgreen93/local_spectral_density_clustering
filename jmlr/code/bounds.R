#----------------------------------------------#
# Theoretical bounds
#----------------------------------------------#

ncut_bound <- function(r,candidate_cluster,distribution){
  # Normalized cut bound (Prop 12.) applied to 
  # mixture of uniforms with epsilon background noise over domain.
  
  # Parameters
  d <- attributes(candidate_cluster)[["d"]]
  sigma <- attributes(candidate_cluster)[["sigma"]]
  nu_cluster <- attributes(candidate_cluster)[["measure"]]
  epsilon <- distribution$epsilon
  nu_domain <- attributes(distribution$domain)[["measure"]]
  n_mixtures <- distribution$n_mixtures
  lambda <- (1 - epsilon)/(n_mixtures * nu_cluster)
  gamma <- 0
  theta <- lambda - epsilon/(nu_domain - 2*nu_cluster)
  
  # TODO: enforce!
  # stopifnot(r <= sigma/(4*d))
  
  # Bound
  (1 - d*r/sigma)^{-2} * d * r/sigma * lambda * (lambda - theta * r^{gamma}/(gamma + 1))/lambda^2
}

conductance_bound <- function(r,candidate_cluster,distribution){
  # Conductance bound (Prop 13.) applied to 
  # mixture of uniforms with epsilon background noise over domain.
  
  # Parameters
  d <- attributes(candidate_cluster)[["d"]]
  sigma <- attributes(candidate_cluster)[["sigma"]]
  L <- attributes(candidate_cluster)[["L"]]
  rho <- attributes(candidate_cluster)[["rho"]]
  
  # TODO: enforce necessary constraints!
  
  # Bound
  sqrt(2*pi)/36 * (1 - r/(4*rho*L)) * (1 - r/sigma * sqrt((d + 2)/(2 * pi)) )^2 * 
    r/(rho*L*sqrt(d + 2))
}

local_spread_bound <- function(r, candidate_cluster, distribution){
  # Local spread bound (Lemma 11.) applied to 
  # mixture of uniforms with epsilon background noise over domain.
  
  # Parameters
  d <- attributes(candidate_cluster)[["d"]]
  sigma <- attributes(candidate_cluster)[["sigma"]]
  L <- attributes(candidate_cluster)[["L"]]
  rho <- attributes(candidate_cluster)[["rho"]]
  
  # Bound
  1/4 * (2*r/rho)^d * (1 - r/sigma * sqrt((d + 2)/(2*pi)))
}