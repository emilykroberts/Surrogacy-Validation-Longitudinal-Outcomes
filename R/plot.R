
#' Plot
#'
#' @description Plot ...
#'
#' @param samp samp ...
#'
#' @return a plot
#'
#' @examples
#'
#' set.seed(323);
#' array_id2 = array_id = 1
#' sim = 2
#'
#' # simulation parameters
#' n = 300; SIM = 1000 # burnin, sample size, and number of MCMC iterations
#' accept = 0 # counter for acceptance ratio
#' condindfit = T # fit conditional independence
#' prior = 'beta' # or unif
#' misspecify = F
#' r = 3
#'
#' # generating values
#' alpha1 = 2; psi2 = 0
#' beta0 = 2; omega1 = 0
#' beta1 = 3.1; omega2 = 0
#' eps1 = 0.5; eps2 = 0.5; eps3 = 0.5; tau4 = 0.5; eta4 = 1
#' theta10 = 0.15; theta11 = 0.7; thetaT = 0.2142857
#' sigmae = 0.3; tau = sigmae^2 # represents ei
#'
#'ST = generate_data(n, theta10 = theta10, theta11 = theta11, thetaT = thetaT,
#'  alpha1 = alpha1, psi2 = psi2,
#'  beta0 = beta0, omega1 = omega1,
#'  beta1 = beta1, omega2 = omega2,
#'  eps1 = eps1, eps2 = eps2, eps3 = eps3, tau4 = tau4, eta4 = eta4,
#'  sigmae = 0.3,
#'  misspecify = misspecify, r = r)
#'
#' # test with true values of random effects
#' params = run_mcmc(samp = samp, Xi = rep(0, n),  n = n, SIM = SIM, r = r, condindfit = condindfit, prior = prior,
#'                   r13init = 0.6, r12init = 0.1, r23init = 0.1,
#'                   b0init = uitrue[1:(n/2),2], b1init = uitrue[(n/2 + 1):n,3])
#'
#' # test with uninformative estimates of random effects
#' params = run_mcmc(samp = samp, Xi = rep(0, n),  n = n, SIM = SIM, r = r, condindfit = condindfit, prior = prior,
#'                   r13init = 0.6, r12init = 0.1, r23init = 0.1,
#'                   b0init = rep(0, n/2), b1init = rep(0, n/2))
#' plot_traceplots(params = params, variable = "int")
#' plot_traceplots(params = params, variable = "slope")
#' plot_traceplots(params = params, variable = "r13")
#' plot_traceplots(params = params, variable = "s1")
#' plot_traceplots(params = params, variable = "s2")
#' plot_traceplots(params = params, variable = "s3")
#' plot_traceplots(params = params, variable = "b1")
#'
plot_traceplots = function(params, variable){
  param = params$params


  plot(eval(parse(text = paste0("param$", variable))), ylab = "Parameter Draw",
       xlab = "MCMC Iteration", main = paste("Traceplot of Parameter", variable))

}
