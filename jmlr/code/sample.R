#---------------------------------------------------------------------#
# Functions to sample data.
#---------------------------------------------------------------------#

make_two_2drectangular_sampler <- function(epsilon,sigma,rho,d){
  function(n){
    X <- matrix(nrow = n,ncol = d)
    for(ii in 1:n)
    {
      if(runif(1) <= epsilon)
      {
        X[ii,] <- runif(d, min = -1, max = 1) # Omega
      }
      else if(runif(1) <= .5)
      {
        X[ii,] <- c(runif(1, min = -1/2 - sigma/2, max = -1/2 + sigma/2),
                    runif(1, min = -rho/2, max = rho/2)) # C_1
      } else{
        X[ii,] <- c(runif(1, min = 1/2 - sigma/2, max = 1/2 + sigma/2),
                    runif(1, min = -rho/2, max = rho/2))  # C_2
      }
    }
    return(X)
  }
}