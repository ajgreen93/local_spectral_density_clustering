#---------------------------------------------------------------------#
# Set related manipulations.
#---------------------------------------------------------------------#

make_partition <- function(sets){
  function(x){
    x_in_set <- sapply(sets,FUN = function(set){set(x)})
    stopifnot(sum(x_in_set) == 1) # Otherwise not a valid partition.
    
    return(sets[[which(x_in_set)]])
  }
}
make_two_2drectangular_cluster_sets <- function(sigma,rho){
  cluster_sets <- list(
    set_1 <- function(x)
    {
      (-1/2 - sigma/2 <= x[1] && x[1] < -1/2 + sigma/2) && (-rho/2 <= x[2] && x[2] < rho/2)
    },
    set_2 <- function(x)
    {
      (1/2 - sigma/2 <= x[1] && x[1] < 1/2 + sigma/2) && (-rho/2 <= x[2] && x[2] < rho/2)
    }
  )
  attributes(cluster_sets[[1]]) <- list("d" = d, "sigma" = sigma, "rho" = rho, "L" = 1, "measure" = rho*sigma) 
  attributes(cluster_sets[[2]]) <- list("d" = d, "sigma" = sigma, "rho" = rho, "L" = 1, "measure" = rho*sigma) 
  return(cluster_sets)
}
