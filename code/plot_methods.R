require(RColorBrewer)

clustering_plotter <- function(clustering_data, adjacency_matrix, data, plot_directory)
{
  # unpack necessary parameters
  for(method in clustering_data)
  {
    clustering <- method$cluster
    color_vector <- clustering_to_color(clustering, data)
    
    plot_filename <- paste0(plot_directory, method$method$method_type, '_cluster.pdf')
    data_plotter(data, adjacency_matrix, color_vector = color_vector, 
                 plot_filename = plot_filename,
                 color_palette = density_cluster_palette)
  }
}

clustering_to_color <- function(clustering, data)
{
  cluster_vector <- clustering
  
  # TODO: make seed node a global parameter. Different algorithms should have the same seed.
  seed_node_manifold_location <- c(.05, -.35)
  d <- ncol(data)
  seed_node_location <- intrinsic_to_extrinsic(d,seed_node_manifold_location)
  seed_node <- which.min(apply(data, 1, FUN = function(x){sum((x - seed_node_location)^2)}))
  
  cluster_vector <- factor(cluster_vector)
  cluster_vector <- relevel(cluster_vector, as.character(cluster_vector[seed_node]) )
  cluster_vector <- relevel(cluster_vector,  '0')
  cluster_vector <- as.numeric(cluster_vector) # convert to right numbers
  
  return(cluster_vector)
}
plotter <- function(cluster_data, adjacency_matrix, data, P = NULL, plot_directory = NULL)
{
  # unpack necessary parameters
  n <- nrow(data)
  
  # compute density
  density_data <- density(P, data)
    
  # compute cluster tree
  cluster_tree_data <- cluster_tree(P, data, density_data = density_data)
  
  ### Plot (1): raw data
  plot_filename <- paste0(plot_directory, 'data', '.pdf')
  data_plotter(data, plot_filename = plot_filename)
  
  ### Plot (2): graph
  plot_filename <- paste0(plot_directory, 'graph', '.pdf')
  data_plotter(data, adjacency_matrix, plot_filename = plot_filename)
  
  ### Plot (3): density heatmap
  plot_filename <- paste0(plot_directory, 'density_heatmap', '.pdf')
  data_plotter(data,adjacency_matrix, plot_filename = plot_filename,
               color_vector = log(density_data$Kth_neighbor_radius) )
  
  ### Plot (4): density cluster
  
  # form density cluster colors
  p <- length(P$cluster_probabilities)
  cluster_tree_evaluations <- unlist(apply(cluster_tree_data$empirical_cluster_tree,1,FUN = function(c){oracle_density_evaluation(c,P)}))
  density_cut <- c(which.min(cluster_tree_evaluations)[1])
  names(density_cut) <- c('density_cut')
  cluster_vector <- tree_to_clustering(cluster_tree_data$empirical_cluster_tree, density_cut, P)
  color_vector <- clustering_to_color(cluster_vector,data)
  
  # find normalized cut of cluster
  full_density_cluster <- which(cluster_vector == '1')
  ncut <- normalized_cut(full_density_cluster,adjacency_matrix)
  
  
  plot_filename <- paste0(plot_directory, 'true_density_cluster', '.pdf')
  data_plotter(data, adjacency_matrix, color_vector = color_vector, 
               plot_filename = plot_filename,
               color_palette = recovered_density_cluster_palette,
               ncut = ncut)
  
  # heatmap and clusters for clustering algorithms
  for(method in cluster_data)
  {
    parameters <- method$parameters
    
    if(method$method$method_type == 'ppr')
    {
      seed_node <- strtoi(parameters['seed_node'])
    } else{
      seed_node <- NULL
    }
    
    for(color_variable in c('cluster'))
    {
      if(color_variable == 'similarity')
      {
        color_vector <- method[[color_variable]] 
        ncut <- NULL
      } else if (color_variable == 'cluster'){
        color_vector <- as.factor(ifelse((1:n) %in% method[[color_variable]], 1, 0))
        ncut <- normalized_cut(method[[color_variable]],adjacency_matrix)
      }
      plot_filename <- paste0(plot_directory, method$method$method_type, '_', color_variable, '.pdf')
      data_plotter(data,adjacency_matrix, color_vector,seed_node, plot_filename,ncut = ncut)
    }
  }
  
 
  
  # cluster tree evaluation plot
  # clustertree_plot_filename <- paste0(plot_directory, 'clustertree.pdf')
  # clustertree_evaluation_plotter(cluster_data, data, P,clustertree_plot_filename, cluster_tree_data)
  # 
  # # density cluster recovery data plots
  # density_cluster_recovery_data_plotter(cluster_data, data, adjacency_matrix, cluster_tree_data)
}

data_plotter <- function(data,adjacency_matrix = NULL, color_vector = NULL, seed_node = NULL, 
                         plot_filename = NULL, color_palette = NULL, legend = NULL, ncut = NULL)
{
  # TODO: Add function description here
  
  # necessary parameters
  n <- nrow(data)
  
  # form coordinates for graph edges
  if(!is.null(adjacency_matrix))
  {
    adj_list <- which(adjacency_matrix == 1, arr.ind = T)
    x_1 <- data[adj_list[,1],1:2]
    x_2 <- data[adj_list[,2],1:2]
  }
  
  # demarcate seed
  if(!is.null(seed_node))
  {
    seed_node_symbol <- ifelse((1:n) == seed_node, 3, 20) # TODO:
    seed_node_size <- ifelse((1:n) == seed_node, 2, 1) # TODO:
  } else{
    seed_node_symbol <- 20
    seed_node_size <- 1
  }
  if(!is.null(color_vector))
  {
    if(!is.null(color_palette))
    {
      color_scheme <- color_palette
      coloring <- color_scheme[color_vector]
      if(!is.null(seed_node))
      {
        coloring[seed_node] <- 'black'
      }
    }
    else if(length(unique(color_vector)) < 7)
    {
      ncolors <- nlevels(color_vector)
      color_scheme <- cluster_color_pallete(ncolors)
      coloring <- color_scheme[color_vector]
      coloring[seed_node] <- 'black'
    } else
    {
      coloring <- numeric(n)
      color_scheme <- similarity_color_pallete(10)
      if(!is.null(seed_node))
      {
        coloring[-seed_node] <- color_scheme[as.numeric(cut(color_vector[-seed_node],breaks = 10))]
        coloring[seed_node] <- 'black'
      } else
      {
        coloring <- color_scheme[as.numeric(cut(color_vector,breaks = 10))]
      }
      
    }
  } else
  {
    coloring <- 'black'
  }
  
  # legend material
  if(!is.null(legend))
  {
    legend_label <- legend[color_scheme %in% unique(coloring)]
    legend_color <- color_scheme[color_scheme %in% unique(coloring)]
  }
  
  
  # plot
  if(!is.null(plot_filename))
  {
    pdf(plot_filename, width = par()$fin[1], height = par()$fin[2])
  }
  old.par <- par(mai = c(0.1,0.1,0.1,0.1))
  plot(x = data[,1],
       y = data[,2],
       pch = 20,
       col = coloring,
       axes = F,
       ann = F,
       type= 'n')
  if(!is.null(adjacency_matrix))
  {
    for(i in 1:nrow(x_1))
    {
      segments(x_1[i,1],x_1[i,2],x_2[i,1],x_2[i,2],col = 'grey')
    }
  }
  points(x = data[,1],
         y = data[,2],
         pch = seed_node_symbol,
         cex = seed_node_size,
         col = coloring)
  if(!is.null(legend))
  {
    legend('topright',
           legend = legend_label,
           fill = legend_color)
  }
  if(!is.null(ncut))
  {
    # title(main = paste0('Cluster normalized cut: ', round(ncut,3)))
  }
  par(old.par)
  if(!is.null(plot_filename))
  {
    dev.off()
  }
  
}

density_cluster_recovery_data_plotter <- function(cluster_data, data, adjacency_matrix,cluster_tree_data)
{
  # unpack necessary parameters
  n <- nrow(data)
  
  # find closest density cluster to actual cluster
  for(method in cluster_data)
  {
    # to put on data plot
    parameters <- method$parameters
    seed_node <- strtoi(parameters['seed_node'])
    
    cluster <- method$cluster
    evaluation_data <- density_clustertree_recovery(cluster, data, P, cluster_tree_data)
    
    # Plot 1 -- high density regions
    
    # color by
    color_vector <- rep(1,n)
    color_vector[evaluation_data$closest_density_cluster] <- 2
    color_vector[evaluation_data$other_density_clusters] <- 3
    
    plot_filename <- paste0(plot_directory, method$method$method_type, 'recovered_density_cluster.pdf')
    data_plotter(data,adjacency_matrix, color_vector, seed_node, 
                 plot_filename, recovered_density_cluster_palette, recovered_density_cluster_legend)
    
    # Plot 2 -- cluster recovery
    
    # color by
    # 1. estimated cluster
    # 2. best density cluster
    # 3. grey for all other points
    color_vector <- rep(1,n)
    color_vector[evaluation_data$closest_density_cluster] <- 2
    color_vector[evaluation_data$other_density_clusters] <- 3
    color_vector[cluster] <- 4
    
    
    # data plot
    plot_filename <- paste0(plot_directory, method$method$method_type, 'density_cluster_recovery.pdf')
    data_plotter(data, adjacency_matrix, color_vector, seed_node, 
                 plot_filename, density_cluster_recovery_palette, density_cluster_recovery_legend)
  }
}

clustertree_evaluation_plotter <- function(cluster_data, data, P,plot_filename, cluster_tree_data)
{
  method_names <- sapply(cluster_data,FUN = function(method){method$method$method_type})
  clustertree_recovery <- matrix(ncol = length(cluster_data) + 1, 
                            nrow = 100,
                            dimnames = list(NULL, 
                                            c('density_threshold',method_names))) # TODO: this is just arbitrarily fixed to be the same length as approx
  for(method in cluster_data)
  {
    cluster <- method$cluster
    evaluation_data <- density_clustertree_recovery(cluster, data, P, cluster_tree_data)
    clustertree_recovery[,method$method$method_type] <- evaluation_data$cluster_tree_evaluation[,'cluster_recovery']
  }
  clustertree_recovery[,'density_threshold'] <- evaluation_data$cluster_tree_evaluation[,'density_threshold']
  
  # plot
  colors = brewer.pal(length(method_names),'Set1')
  names(colors) <- method_names
  if(!is.null(plot_filename))
  {
    pdf(plot_filename)
  }
  plot(x = clustertree_recovery[,'density_threshold'],
       y = seq(0,1,length.out = nrow(clustertree_recovery)),
       type = 'n')
  for(method_name in method_names)
  {
    lines(x = clustertree_recovery[,'density_threshold'],
          y = clustertree_recovery[,method_name],
          col = colors[method_name])
  }
  dev.off()
  
}

density_cluster_palette <- c("grey30", "red", 'blue')
recovered_density_cluster_palette <- c("grey30", "orange", 'cyan1')
recovered_density_cluster_legend <- c('no_cluster', 'density_cluster','density_level_set')

density_cluster_recovery_palette <- c("grey30", "blue", "red", 'green4')
density_cluster_recovery_legend <- c('no_cluster', 'density_cluster', 'density_level_set','estimated_cluster')

cluster_color_pallete <- colorRampPalette(
  colors = c("blue", "orange", "red"),
  space = "Lab" # Option used when colors do not represent a quantitative scale
)

similarity_color_pallete <- colorRampPalette(
  colors = c('blue','red')
)

### Function for saving pdfs without margins
savepdf <- function(file, width=16, height=10)
{
  fname <- file
  pdf(fname, width=width/2.54, height=height/2.54,
      pointsize=10)
  par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.1,1.1))
}


