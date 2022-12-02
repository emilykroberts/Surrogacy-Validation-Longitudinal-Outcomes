#' traceplots
#'
#' @description creates traceplots
#'
#' @param params matrix of parameters from mcmc
#' @param variable variable of which to plot
#'
#' @return traceplots
#'
#' @examples
#' example(plot_traceplots(params = params, variable = "int"))

#'
plot_traceplots = function(params, variable){
  param = params$params

  plot(eval(parse(text = paste0("param$", variable))), ylab = "Parameter Draw",
       xlab = "MCMC Iteration", main = paste("Traceplot of Parameter", variable))

}
