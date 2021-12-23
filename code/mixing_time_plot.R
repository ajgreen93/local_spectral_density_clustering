#----------------------#
# script to generate plot for mixing time data
#----------------------#

library(sfsmisc)

# load data
source_path <- 'C:/Users/alden/OneDrive/Documents/Statistics/local_spectral_density_clustering/data/'
data_directory <- 'mixing_time_data/sigma_data/'
data_path <- paste0(source_path, data_directory)
for(file in list.files(data_path))
{
  load(paste0(data_path,file))
}

# setup plot file
plot_path <- paste0(source_path, 'example1plots/')
plot_filename <- paste0(plot_path,'sigma_mixing_time_plot.pdf')
pdf(plot_filename, 6, 6)

# compute necessary quantities for plot
theoretical_data <- theoretical_mixing_time
empirical_data <- empirical_mixing_time

plot_empirical_data <- t(apply(empirical_data,c(1,2),FUN = function(x){median(x,na.rm = T)}))
plot_theoretical_data <- matrix(0, nrow = nrow(theoretical_data), ncol = ncol(theoretical_data))
for(ii in 1:nrow(theoretical_data))
{
  plot_theoretical_data[ii,] <- theoretical_data[ii,] * 1.2*(max(plot_empirical_data[ii,]/theoretical_data[ii,]))
}

sigma <- configs$model$sigma
x_var_range <- configs$model$sigma_range + 2 * sigma
# plot_df <- data.frame(x = log(D_range), 
#                       y_emp = log(med_empirical_mixing_time),
#                       y_the = log(t(theoretical_mixing_time))
#                       )
col <- c('green4','mediumpurple4')
# col <- c('chartreuse','deeppink')

# plot specifications
line.lwd <- 2
y.min <- min(c(plot_empirical_data, plot_theoretical_data))
y.max <- max(c(plot_empirical_data, 
               plot_theoretical_data))
ylim <- c(y.min, y.max)
old.par <- par(no.readonly = TRUE)
x.at <- x_var_range

# shadow plot
plot(x_var_range, plot_empirical_data[1,], log = 'xy',
     xlab = 'sigma', ylab = 'Mixing Time', ylim = ylim, axes = FALSE, frame.plot = TRUE, type = 'n',
     cex.lab = 1.7)
x_axis <- axis(1, cex.axis = 1.7)
rug_x <- numeric(length(x_axis) - 1)*4
for(ii in 2:(length(x_axis)))
{
  for(jj in 1:4)
  {
    kk <- 4*(ii - 2) + jj
    rug_x[kk] <- x_axis[ii] - (x_axis[ii] - x_axis[ii - 1])/(2^jj)
  }
}
rug(rug_x, ticksize = -.015, side = 1)
y_axis <- axis(2,cex.axis = 1.7)
rug_y <- numeric(length(y_axis) - 1)*4
for(ii in 2:(length(x_axis)))
{
  for(jj in 1:4)
  {
    kk <- 4*(ii - 2) + jj
    rug_y[kk] <- y_axis[ii] - (y_axis[ii] - y_axis[ii - 1])/(2^jj)
  }
}
rug(rug_y, ticksize = -.015, side = 2)


# plot each empirical line
old.par <- par(lwd = line.lwd, no.readonly = TRUE)
empirical_lty <- 2
empirical_pch <- 2

lines(x_var_range, plot_empirical_data[1,], lty = empirical_lty,col = col[1])
lines(x_var_range, plot_empirical_data[2,], lty = empirical_lty,col = col[2])
points(x_var_range, plot_empirical_data[1,], pch = empirical_pch,col = col[1])
points(x_var_range, plot_empirical_data[2,], pch = empirical_pch,col = col[2])

# plot each theoretical line
theoretical_lty <- 1
theoretical_pch <- 4
lines(x_var_range, plot_theoretical_data[1,], lty = theoretical_lty, col = col[1])
lines(x_var_range, plot_theoretical_data[2,], lty = theoretical_lty, col = col[2])
points(x_var_range, plot_theoretical_data[1,], pch = theoretical_pch, col = col[1])
points(x_var_range, plot_theoretical_data[2,], pch = theoretical_pch, col = col[2])

# back to the old line width
par(old.par)

# add grid
grid(equilogs = FALSE, lwd = 1)

# add legend
leg.txt <- c("2d empirical", "2d theoretical", "3d empirical", "3d theoretical")
leg.col <- c(col[1],col[1],col[2],col[2])
leg.lty <- c(2,1,2,1)
leg.pch <- c(2,4,2,4)
legend("topright", leg.txt, col = leg.col, lty = leg.lty, pch = leg.pch, inset = .02, cex = 1.2)

# finish plot
dev.off()