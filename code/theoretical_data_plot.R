# theoretical_data_plot.R is the main script
# for plotting the data, complete with highlighting of the density cluster and expanded set




# Dependencies
require(MASS)
require(ggplot2)
require(igraph)
require(RANN)
require(Matrix)
require(reshape2)

# Source methods
source_path <- 'C:/Users/alden/OneDrive/Documents/Statistics/local_spectral_density_clustering/'
code_source_path <- paste(source_path, "code/", sep = '')
data_path <- paste0(source_path, "data/")
script_names <- c('configs/normalized_cut_config.R','theoretical_methods.R','sample.R',
                  'geometric_methods.R','numeric_methods.R','graph.R')
for(script_name in script_names)
{
  source(paste(code_source_path, script_name, sep = ''))
}

# compute theoretical and sample normalized cut

# unpack necessary parameters
distribution_class <- model$distribution_class
sigma_range <- model$sigma_range
gamma <- model$gamma
d_range <- model$d_range
lambda_sigma <- model$lambda_sigma
lambda <- model$lambda

#---------------------#
# Plot 1: sigma = .466, d = 2
#---------------------#

sigma <- sigma_range[5]
d <- d_range[1]
P <- list(
  distribution_class = distribution_class,
  sigma = sigma,
  gamma = gamma,
  D = D,
  d = d,
  lambda_sigma = lambda_sigma,
  lambda = lambda
)

file_path <- paste0(data_path, "sample1.pdf")

plot_sample(P,n,file_path)

#---------------------#
# Plot 1: sigma = 3.2, d = 2
#---------------------#

sigma <- sigma_range[10]
d <- d_range[1]
P <- list(
  distribution_class = distribution_class,
  sigma = sigma,
  gamma = gamma,
  D = D,
  d = d,
  lambda_sigma = lambda_sigma,
  lambda = lambda
)

file_path <- paste0(data_path, "sample2.pdf")

plot_sample(P,n,file_path)


### the plotting function

plot_sample <- function(P,n, filename = NULL)
{
  if(!is.null(filename))
  {
    pdf(filename)
  }
  
  ## Sample data
  X <- distribution_sample(P,n)
  
  ## Compute density cluster and expanded set
  
  D <- P$D
  sigma <- P$sigma
  
  C <- I_Csig(X,D,sigma,0)
  Csig <- I_Csig(X,D,sigma,sigma)
  
  cluster_color <- ifelse((1:n) %in% C, 'red', 
                          ifelse((1:n) %in% Csig,'orange','blue'))
  
  plot(X, col = alpha(cluster_color,.16), asp = 1, pch = 20,
       xlab = 'x1', ylab = 'x2', cex.lab = 1.7, cex.axis = 1.7)
  
  if(!is.null(filename))
  {
    dev.off()
  }
}

