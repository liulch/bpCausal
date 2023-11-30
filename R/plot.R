plot_function <- function(x, y, treatment = 0, 
                          dataset = NULL, model, shadow = TRUE) {
  if (is.null(x)) {
    warning("x is not the correct. x can be something like x <- c(-19:10)")
  }
  
  if (is.null(y)) {
    warning("y is not the correct. y can be something like y <- apply(matrix(simdata[which(simdata$treat==1),'eff'], 30, 5), 1, mean)")
  }
  
  if (!is.numeric(y)) {
    warning("y1 is not of type numeric.")
  }
  
  if(!is.numeric(treatnent)) {
    warning("The value of treatment should be a number.")
  }
  
  x_range <- range(x)
  if(x > x_range[2] || x < x_range[1]) {
    warning("The value of treatment must be in the range of x.")
  }
  
  if (!class(model) == "bpCausal"){
    warning("The model is not a bpCausal model.")
  }
  
  if (shadow != TRUE || shadow != FALSE) {
    warning("shadow should either be TRUE or FALSE.")
  }
  
 
  color <-rgb( 0.06, 0,0.1,0.1)
  est_ATT_range <- range(eout1$est.eff$estimated_ATT)
  plot_x_lim <- c(x_range[1],x_range[2])
  plot_y_lim <- c(ceiling(est_ATT_range[1] - 1), ceiling(est_ATT_range[2] + 1))
  eout1 <- effSummary(model)
  estimated_ATT <- eout1$est.eff$estimated_ATT
  estimated_ATT_ci_l <- eout1$est.eff$estimated_ATT_ci_l
  estimated_ATT_ci_u <- eout1$est.eff$estimated_ATT_ci_u
  
  plot(x, y, type = "l", col = "red", ylim = plot_y_lim, 
       xlab = "Time", ylab = "Estimated ATT", cex.lab = 1.5)
  abline(v = treatment, lty = 3, col = "grey")

  lines(x1, estimated_ATT)
  lines(x1, estimated_ATT_ci_l, lty = 2)
  lines(x1, estimated_ATT_ci_u, lty = 2)
  
  if (shadow == TRUE) {
    polygon(c(x1, rev(x1)), c(estimated_ATT_ci_l, rev(estimated_ATT_ci_u)), 
            col = color, border = NA)
  }
 
}