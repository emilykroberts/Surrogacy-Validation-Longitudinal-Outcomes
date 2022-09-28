
#' Generate data
#'
#' @description Generate data ...
#'
#' @param n number of observations to simulate
#'
#' @return a list of samples
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
#' samp = ST[[1]]
#' uitrue = ST[[2]]
#'
generate_data = function(n, theta10, theta11, thetaT, alpha1, psi2,
                         beta0, omega1,
                         beta1, omega2,
                         eps1, eps2, eps3, tau4, eta4,
                         sigmae, misspecify, r){

  tau = sigmae^2
  R = matrix(rep(1,3*3),3,3); R[1,2] = R[2,1] = theta10; R[1,3] = R[3,1] =theta11; R[2,3] = R[3,2] = thetaT
  S = diag(c(eps1, eps2, eps3));

  ### generate longitudinal T0 and T1 = mu + ui + eij #ui shared across outcomes per person
  biR = R

  allsamp = mvrnorm(n, mu = c(alpha1, rep(c(beta0, beta1), r)) , diag(c(rep(tau, 2*r + 1))))
  # update to make more general

  allsamp[,1] = alpha1

  # allow for misspecifying MVN distribution of random effects
  uitruesamp = uitrue = mvrnorm(n, c(0, 0, 0), S %*% biR %*% S)

  if(misspecify){
    uitrue = rmvgamma(n, shape = 0.5, rate = 0.5, corr = S %*% biR %*% S)
    uitruesamp = uitrue = uitrue - colMeans(uitrue)
  }

  allsamp[, 1] = allsamp[, 1] + uitrue[,1]

  allsamp[, 1 + seq(1, 2*r, 2)] = allsamp[, 1 + seq(1, 2*r, 2)] + uitrue[,2]
  allsamp[, 2 + seq(1, 2*r, 2)] = allsamp[, 2 + seq(1, 2*r, 2)] + uitrue[,3]

  allsamp = cbind(allsamp, rep(0, n))
  allsamp[(1:n/2), seq(1, 2*r + 2, 2)] = NA
  allsamp[(n/2+1):n, 1 + seq(1, 2*r , 2)] = NA # update to make more general


  allsamp_list = list(allsamp, uitrue)
  return(allsamp_list)

}
