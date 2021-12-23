library(fields)

## input parameters
pdf <- TRUE
filename <- "lower_bound_1"
rho <- 1
sigma <- .2
epsilon <- .2
n_eval <- 200^2
x_range <- c(-.6,.6)
y_range <- c(-.6,.6)

density_maker <- function(rho,sigma,epsilon)
{
  density <- function(x)
  {
    if(x[1] < -3*sigma/2 || x[1] > 3*sigma/2 || x[2] < -rho/2 || x[2] > rho/2)
    {
      return(0)
    }
    else if(x[1] <= -sigma/2 || x[1] >= sigma/2)
    {
      return( (1-epsilon)/(2*sigma*rho) )
    } else{
      return(epsilon/(sigma*rho))
    }
  }
}

density <- density_maker(rho,sigma,epsilon)
eval_points <- matrix(0,ncol = 2,nrow = n_eval)
for(ii in 1:sqrt(n_eval))
{
  x <- x_range[1] + ii/sqrt(n_eval)*(x_range[2] - x_range[1])
  for(jj in 1:sqrt(n_eval))
  {
    kk <- (ii - 1)*sqrt(n_eval) + jj
    y <- y_range[1] + jj/sqrt(n_eval)*(y_range[2] - y_range[1])
    eval_points[kk,] <- c(x,y)
  }
}
p <- apply(eval_points,1,FUN = density)

if(pdf)
{
  source_path <- 'C:/Users/alden/OneDrive/Documents/Statistics/local_spectral_density_clustering/'
  data_path <- paste0(source_path, "data/")
  file_path <- paste0(data_path,filename,".pdf")
  pdf(file_path)
}

# create colors
n_colors <- 100
colfunc <- colorRampPalette(c("blue", "red"))
colors <- colfunc(n_colors)
color_range <- c(.25,2.25)
# breaks <- seq(min(color_range),max(color_range),length.out = n_colors)
# p_colors <- colors[cut(p,breaks,include.lowest = T)]
# plot(eval_points, pch = 20, col = p_colors,
#       xlab = 'x1', ylab = 'x2', xlim = c(-.6,.6), ylim = c(-.6,.6), cex.lab = 1.2, cex.axis = 1.2)

breaks <- seq(min(color_range),max(color_range),length.out = n_colors + 1)
plot_data <- matrix(p,nrow = sqrt(n_eval),ncol = sqrt(n_eval),byrow = TRUE)
image.plot(plot_data, nlevel = n_colors, col = colfunc(n_colors), breaks = breaks,
           legend.cex = 1.5,axes = F)
axis(2, at = seq(0,1,length.out = 6), labels = seq(y_range[1],y_range[2],length.out = 6),cex.axis = 1.5)
axis(1, at = seq(0,1,length.out = 6), labels = seq(y_range[1],y_range[2],length.out = 6), cex.axis = 1.5)
box()
if(pdf)
{
  dev.off()
}