plot_att <- function(x, y, treatment = 0, eout, shadow = TRUE) {
  
  x_range <- range(x)
  
  color <-rgb( 0.06, 0,0.1,0.1)
  est_ATT_range <- range(eout1$est.eff$estimated_ATT)
  plot_x_lim <- c(x_range[1],x_range[2])
  plot_y_lim <- c(ceiling(est_ATT_range[1] - 1), ceiling(est_ATT_range[2] + 1))
  estimated_ATT <- eout1$est.eff$estimated_ATT
  estimated_ATT_ci_l <- eout1$est.eff$estimated_ATT_ci_l
  estimated_ATT_ci_u <- eout1$est.eff$estimated_ATT_ci_u
  
  plot(x, y, type = "l", col = "red", ylim = plot_y_lim, 
       xlab = "Time", ylab = "Estimated ATT", cex.lab = 1.5)
  abline(v = treatment, lty = 3, col = "grey")
  
  lines(x1, estimated_ATT, col = "gray33")
  lines(x1, estimated_ATT_ci_l, lty = 2, col = "darkgrey")
  lines(x1, estimated_ATT_ci_u, lty = 2, col = "darkgrey")
  
  if (shadow == TRUE) {
    polygon(c(x1, rev(x1)), c(estimated_ATT_ci_l, rev(estimated_ATT_ci_u)), 
            col = color, border = NA)
  }
  
}
