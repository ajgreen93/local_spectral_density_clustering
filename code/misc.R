plot_df <- data.frame(x_1 = g_enl[,1], x_2 = g_enl[,2], col = q_dens)

ggplot(plot_df, aes(x = x_1, y = x_2, color = col)) + geom_point()

plot(X, col = alpha(ifelse((1:n) %in% Csig, 'red', 'blue'),.04), asp = 1, pch = 20,
     xlab = 'x1', ylab = 'x2')

plot(X, col = ifelse((1:n) %in% Csig[rowSums(Csig_adjacency_matrix) == 1], 'red', 'blue'), asp = 1)

C_box <- which(X[,1] <= sigma & X[,1] >= 0 & 
               X[,2] <= D & X[,2] >= 0)


plot(X, col = ifelse((1:n) %in% C_box, 'red', 'blue'), asp = 1)


plot(X, col = alpha(ifelse((1:n) %in% Csig, 'red', 'blue'),.04), asp = 1, pch = 20,
     xlab = 'x1', ylab = 'x2')

