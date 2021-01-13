#----------------------------------------------------#
# Plot results from validate theory pipeline
#----------------------------------------------------#

# What is to be plotted? # (NOTE: this script is not robust to changing these; they are global variables.)
vectors <- c("samples","normalized_ppr","sweep_cut")
quantities <- c("normalized_cut","conductance","local_spread","normalized_volume_ssd")
n_vectors <- 3
n_quantities <- 4
parameter <- "diameter"

# Please specify the name of the directory containing the saved data here.
dir_name <- "20210111200827"

# Would you like to save plots to file?
save_plots <- TRUE

# Which distributions (by index) would you like to plot
plot_indx <- 1:10

# Structure for storing data to plot
list_of_plot_data <- vector(mode = "list", length = n_quantities)
names(list_of_plot_data) <- quantities

# Load results
save_directory <- paste0("data/",dir_name)
for(file_name in list.files(save_directory, pattern = "*.R"))
{
  load(paste0(save_directory,"/",file_name))
}
empirical_data <- list(ncut = empirical_ncuts,
                       conductance = empirical_conductances,
                       local_spread = empirical_local_spreads,
                       normalized_volume_ssd = empirical_normalized_volume_ssd)
population_bounds <- list(ncut = population_ncuts,
                       conductance = population_conductances,
                       local_spread = population_local_spreads,
                       normalized_volume_ssd = population_condition_numbers)

# Data for plotting
distributions <- configs$distributions
seed_location <- configs$seed_location
rhos <- sapply(distributions,FUN = function(distribution){
  attributes(distribution$partition(seed_location))[["rho"]]
  }
)
for(ii in 1:n_quantities)
{
  avg_empirical <- rowMeans(empirical_data[[ii]])[plot_indx]
  list_of_plot_data[[ii]] <- list(empirical = avg_empirical,
                                  population = population_bounds[[ii]][plot_indx])
} 

# Generic plotting functions.
generic_2ddata_plot <- function(X,colors,pch){
  plot(X[,1],X[,2],col = colors, pch = 16, xlab = "x1", ylab = "x2")
}

generic_quantity_plot <- function(x, list_of_ys, log = "xy", xlab = "diameter", ylab = "y",
                      colors, pchs, lwd = 2, cex = cex){
  ylims <- c(.5 * min(unlist(list_of_ys)), 2 * max(unlist(list_of_ys)))
  xlims <- c(min(x),max(x))
  plot(x = x,y = list_of_ys[[1]], ylim = ylims, type = "n",log = log, xlab = xlab,ylab = ylab,
       cex.lab = cex, cex.axis = cex)
  grid(lwd = 3*lwd)
  for(ii in 1:length(list_of_ys))
  {
    y <- list_of_ys[[ii]]
    lines(x = x,
          y = y,
          col = colors[ii],
          lwd = lwd)
    points(x = x,
           y = y,
           col = colors[ii],
           pch = pchs[ii],
           cex = cex)
  }
}

color.gradient <- function(x, colors=c("blue","red"), colsteps=100, transform = function(x){x}) {
  return( colorRampPalette(colors) (colsteps) [ findInterval(transform(x), seq(quantile(transform(x),.01),quantile(transform(x),.995), length.out=colsteps)) ] )
}

# Generate some plots
if(save_plots){
  plot_dir <- file.path("data",dir_name,"plots")
  if(!dir.exists(plot_dir)) dir.create(plot_dir)
}

# Plot the samples, the normalized (by degree) PPR vector, and resulting sweep cut
ppr_color_transform <- sqrt
for(ii in plot_indx)
{
  X <- Xs[[ii]]
  p <- ps[ii,]
  estimated_cluster <- estimated_clusters[ii,]
  
  # TODO: put this is main pipeline
  library(magrittr)
  library(Matrix)
  library(RANN)
  library(reshape2)
  source("graph.R")
  r <- configs$r
  G <- neighborhood_graph(X,r)
  p_norm <- ifelse(is.na(p/degree(G)),0,p/degree(G))
  
  # Colors
  n <- nrow(X)
  colors_for_each_plot <- vector(mode = "list", length = 3)
  names(colors_for_each_plot) <- vectors
  
  colors_for_each_plot[["samples"]] <- rep("black",n)          
  colors_for_each_plot[["normalized_ppr"]] <- color.gradient(p_norm, transform = ppr_color_transform)
  colors_for_each_plot[["sweep_cut"]] <- ifelse(estimated_cluster,"red","blue")
  
  for(plot_name in vectors){
    if(save_plots)
    {
      file_name <- paste0(plot_name,ii,".pdf")
      plot_path <- file.path(plot_dir,file_name)
      pdf(plot_path, 8, 8)
    }
    generic_2ddata_plot(X,colors_for_each_plot[[plot_name]])
    if(save_plots){dev.off()}
  }
}

# Plot the quantities we use to bound volume of SSD
colors <- c("red","blue")
pchs <- c(15,16)
log <- "xy"
if(save_plots){
  lwd <- 2
  cex <- 2
} else{
  lwd <- 1
  cex <- 1
}
for(ii in 1:n_quantities)
{
  if(save_plots)
  {
    file_name <- paste0(quantities[ii],".pdf")
    plot_path <- file.path(plot_dir,file_name)
    pdf(plot_path, 8, 8)
  }
  generic_quantity_plot(rhos, list_of_plot_data[[ii]], xlab = parameter, ylab = names(list_of_plot_data)[ii],
               colors = colors, pchs = pchs, lwd = lwd, cex = cex, log = log)
  if(save_plots)
  {
    dev.off()
  }
}


