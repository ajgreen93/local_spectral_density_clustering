load(file = paste0(save_path, 'metric_evaluation_table.R'))
metric_evaluation_table_isotropic <- metric_evaluation_table

load(file = paste0(save_path, 'metric_evaluation_table_2.R'))
metric_evaluation_table_anisotropic <- metric_evaluation_table

symmetric_set_difference <- matrix(c(rowMeans(metric_evaluation_table_isotropic),
                                     rowMeans(metric_evaluation_table_anisotropic)), ncol = 2)
cluster_center_distance <- P_sequence$mu_distance_scale * .5

plot_df = data.frame(X = cluster_center_distance, 
                     Y_1 = symmetric_set_difference[,1],
                     Y_2 = symmetric_set_difference[,2])
plot_df <- melt(plot_df, id = 'X')
example3_plot1 <- ggplot(data = plot_df, aes(x = X, y = value, colour = variable, shape = variable)) + 
  geom_line(size = 1.2) +
  scale_colour_manual(name = 'Mixture model', labels = c('GMM-isotropic', 'GMM-anisotropic'),
                      values = c('orangered1', 'purple')) +
  geom_point(size = 5) +
  scale_shape_manual(name = 'Mixture model', labels = c('GMM-isotropic', 'GMM-anisotropic'),
                     values = c(18, 20)) +
  xlab('Distance between component centers') + ylab('symmetric set difference') + 
  theme_classic() +
  theme(legend.justification=c(1,1), legend.position=c(1,1), 
        legend.background = element_rect(linetype = 1, size = 0.5, colour = 1),
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16))