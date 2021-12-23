### load each of the four graphs
load(file = paste0(save_path, 'example2_plot1.R'))
load(file = paste0(save_path, 'example2_plot2.R'))
load(file = paste0(save_path, 'example2_plot3.R'))
load(file = paste0(save_path, 'example2_plot4.R'))
### arrange them on a panel

pdf(paste0(save_path,'graphic_1.pdf'))
grid.arrange(example2_plot1,
             example2_plot4,
             example2_plot2,
             example2_plot3,
             nrow = 2)
dev.off()
