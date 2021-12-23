#----------------------------------------------------#
# Plot results from validate theory pipeline
# 
# Plotting options are 
# 
# -- vectors <- c("samples","normalized_ppr","sweep_cut")
# -- quantities <- c("normalized_cut","conductance","local_spread","normalized_volume_ssd")
# -- parameter <- c()
#----------------------------------------------------#

# What is to be plotted?
vectors <- "samples"
quantities <- c("normalized_cut", "conductance","local_spread")
parameter <- "sigma"

# Would you like to save plots to file?
save_plots <- TRUE

# Which distributions (by index) would you like to plot
plot_indx <- 1:10

# Structure for storing data to plot
n_quantities <- length(quantities)
list_of_plot_data <- vector(mode = "list", length = n_quantities)
names(list_of_plot_data) <- quantities

# Load results
save_file <- paste0("2d_uniform_rectangles_",parameter,".rda")
load(file.path("data",save_file))
empirical_data <- list(normalized_cut = empirical_ncuts,
                       conductance = empirical_conductances,
                       local_spread = empirical_local_spreads,
                       normalized_volume_ssd = empirical_normalized_volume_ssd)
population_bounds <- list(normalized_cut = population_ncuts,
                       conductance = population_conductances,
                       local_spread = population_local_spreads,
                       normalized_volume_ssd = population_condition_numbers)

# Take only the asked-for quantities.
empirical_data <- empirical_data[quantities]
population_bounds <- population_bounds[quantities]

# Data for plotting
parameter_values <- get(parameter)
if(parameter == "epsilon"){
  parameter <- "density_ratio"
  parameter_values <- 2*rho*sigma/(4 - 2*rho*sigma) * epsilon/(1 - epsilon)
}

for(ii in 1:n_quantities)
{
  avg_empirical <- rowMeans(empirical_data[[ii]])[plot_indx]
  list_of_plot_data[[ii]] <- list(empirical = avg_empirical,
                                  population = population_bounds[[ii]][plot_indx])
} 

# Generic plotting functions.
generic_2ddata_plot <- function(X,colors,subsample = FALSE,cex = 2){
  if(subsample){
    n <- nrow(X)
    indx <- sample(1:n,size = 8000)
    X <- X[indx,]
    colors <- colors[indx]
  }
  plot(X[,1],X[,2],col = colors, pch = 16, xlab = "x1", ylab = "x2",
       cex.lab = cex, cex.axis = cex)
}

generic_quantity_plot <- function(x, list_of_ys, log = "xy", xlab = "diameter", ylab = "y",
                      colors, pchs, lwd = 2, cex = cex){
  ylims <- c(.5 * min(unlist(list_of_ys)), 2 * max(unlist(list_of_ys)))
  xlims <- c(min(x),max(x))
  xlab <- gsub("_", " ", xlab)
  ylab <- gsub("_", " ", ylab)
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
  plot_dir <- file.path("plots",parameter)
  if(!dir.exists(plot_dir)) dir.create(plot_dir,recursive = TRUE)
}

# Plot the samples, the normalized (by degree) PPR vector, and resulting sweep cut
subsample <- TRUE
plot_vectors <- if(is.null(vectors)) FALSE else TRUE
if(plot_vectors){
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
      generic_2ddata_plot(X,colors_for_each_plot[[plot_name]],subsample)
      if(save_plots){dev.off()}
    }
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

plot_quantities <- if(is.null(quantities)) FALSE else TRUE
if(plot_quantities){
  for(ii in 1:n_quantities)
  {
    if(save_plots)
    {
      file_name <- paste0(quantities[ii],".pdf")
      plot_path <- file.path(plot_dir,file_name)
      pdf(plot_path, 8, 8)
    }
    generic_quantity_plot(parameter_values, list_of_plot_data[[ii]], xlab = parameter, ylab = names(list_of_plot_data)[ii],
                          colors = colors, pchs = pchs, lwd = lwd, cex = cex, log = log)
    if(save_plots)
    {
      dev.off()
    }
  }
}


