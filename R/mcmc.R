
#' Run MCMC
#'
#' @description Run MCMC ...
#'
#' @param samp data
#' @param Xi covariates
#' @param n sample size
#' @param SIM number of mcmc iterations
#' @param r number of repeated measures
#' @param condindfit logical value to assume conditional independence
#' @param prior either uniform or beta prior on theta
#' @param r12init initial value of theta10
#' @param r13init initial value of theta11
#' @param r23init initial value of thetaT
#'
#' @return a list of samples
#'
#' @examples
#' example(run_mcmc(samp = samp, Xi = rep(0, n),  n = n, SIM = SIM, r = r, condindfit = condindfit, 
#' prior = prior, r13init = 0.6, r12init = 0.1, r23init = 0.1,
#' b0init = rep(0, n/2), b1init = rep(0, n/2)))

run_mcmc = function(samp, Xi, n, SIM, r, condindfit, prior, r13init, r12init, r23init,
                    b0init, b1init){
  burnin = SIM * .3
  trt = c(rep(0,n)); trt[(n-(n/2-1)):(n)] = 1

  # set up matrices to save results
  holdmu = matrix(rep(0,3*SIM),3,SIM); holdS = holdR = array(rep(0,3*3*SIM),dim = c(3,3,SIM));
  holdR[,,1] = diag(c(1,1,1)); holdS[,,1] = diag(c(1,1,1))
  holdpsi1 = holdpsi2 = holdomega1 = holdomega2 = (rep(0,1*SIM))
  holdalpha1 = holdbeta0 = holdbeta1 = (rep(0,1*SIM));
  int = slope = array(0,c((SIM),1))
  holdB0 = holdB1 = holdSigE = array(0,c((SIM),1))

  bivec = c(b0init, b1init) # can set different starting values of random effects
  ST = cbind(samp[, 1:(2*r + 1)], 0, 0)

  # initial values
  holdR[3,1,1] = holdR[1,3,1] = r13init
  holdR[2,1,1] = holdR[1,2,1] = r12init
  holdR[3,2,1] = holdR[2,3,1] = r23init
  holdB0[1] = holdB1[1] = 1
  holdalpha1[1] = 0; holdbeta0[1] = 0; holdbeta1[1] = 0
  a = b = 0.1

  # start simulation
  sim = 2

  while(sim <= SIM){

    # estimate error term
    tmp5 = c(ST[trt == 1, seq(3, 2*r + 1, 2)] - c(bivec[trt == 1])) - holdbeta1[sim-1]
    tmp4 = c(ST[trt == 0, seq(2, 2*r + 1, 2)] - c(bivec[trt == 0]))  - holdbeta0[sim-1]

    B0 = B1 = rinvgamma(1, shape = a + length(c(tmp4, tmp5))/2, scale = (sum((tmp4 + tmp5)^2)/2 + b))

    tmp5 = c(ST[trt == 1,3] - c(bivec[trt == 1])) - holdbeta1[sim-1]
    tmp4 = c(ST[trt == 0,2] - c(bivec[trt == 0])) - holdbeta0[sim-1]

    # B0 = B1 = rinvgamma(1, shape = a + length(c(tmp4, tmp5))/2, scale = (sum((tmp4 + tmp5)^2)/2 + b))

    holdB0[sim] = B0; holdB1[sim] = B1
    holdSigE[sim] = sqrt(B0)
    Sigma = diag(c(rep(holdSigE[sim]^2, 3)))

    # estimate random effects
    xbeta_all2 = ST
    xbeta_all2[,1] = holdalpha1[sim-1]
    xbeta_all2[,  seq(2, 2*(r), 2)] = holdbeta0[sim-1] # check
    xbeta_all2[,1 + seq(2, 2*(r), 2)] = holdbeta1[sim-1]

    ymisb01 = ymisb11 = ymiss = rep(NA, n)

    # # update to make more general
    S2 = cbind(holdS[,,sim-1] %*% holdR[,,sim-1] %*% holdS[,,sim-1],   matrix(rep(0, 2*r*3), nrow = 3))
    S2 = rbind(S2, matrix(rep(0, r*2*ncol(S2)), ncol = ncol(S2)))
    S2[1, 1 + 1:r] = (holdS[,,sim-1] %*% holdR[,,sim-1] %*% holdS[,,sim-1])[1,2] # cov(S1, T0) 2:4 for r = 3
    S2[1, r + 1 + 1:r] = (holdS[,,sim-1] %*% holdR[,,sim-1] %*% holdS[,,sim-1])[1,3] # cov(S1, T1) 5:7 for r = 4
    S2[1 + 1:r, 1 + 1:r] = (holdS[,,sim-1] %*% holdR[,,sim-1] %*% holdS[,,sim-1])[2,2]
    S2[r + 1 + 1:r, r + 1 + 1:r] = (holdS[,,sim-1] %*% holdR[,,sim-1] %*% holdS[,,sim-1])[3,3]
    S2[1 + 1:r, r + 1 + 1:r] = (holdS[,,sim-1] %*% holdR[,,sim-1] %*% holdS[,,sim-1])[2,3] # cov(T0, T1) = cov(b0, b1)

    for(i in 1:r){
      S2[i + 1, i + 1] = (holdS[,,sim-1] %*% holdR[,,sim-1] %*% holdS[,,sim-1])[2,2] + holdSigE[sim]^2
      S2[i + 1,ncol(S2) - 1] = (holdS[,,sim-1] %*% holdR[,,sim-1] %*% holdS[,,sim-1])[2,2]
      S2[i + 1, ncol(S2)] = (holdS[,,sim-1] %*% holdR[,,sim-1] %*% holdS[,,sim-1])[2,3]

    }
    for(i in 1:r){
      S2[i + 1 + r,i + 1 + r] = (holdS[,,sim-1] %*% holdR[,,sim-1] %*% holdS[,,sim-1])[3,3] + holdSigE[sim]^2
      S2[i + 1 + r,(ncol(S2) - 1)] = (holdS[,,sim-1] %*% holdR[,,sim-1] %*% holdS[,,sim-1])[2,3]
      S2[i + 1 + r,ncol(S2)] = (holdS[,,sim-1] %*% holdR[,,sim-1] %*% holdS[,,sim-1])[3,3]
    }

    S2[1,(ncol(S2) - 1)] = (holdS[,,sim-1] %*% holdR[,,sim-1] %*% holdS[,,sim-1])[1,2]
    S2[1,(ncol(S2))] = (holdS[,,sim-1] %*% holdR[,,sim-1] %*% holdS[,,sim-1])[1,3]

    S2[(ncol(S2) - 1),ncol(S2)] = (holdS[,,sim-1] %*% holdR[,,sim-1] %*% holdS[,,sim-1])[2,3] # b01 b02
    S2[(ncol(S2) - 1),(ncol(S2) - 1)] = (holdS[,,sim-1] %*% holdR[,,sim-1] %*% holdS[,,sim-1])[2,2]
    S2[ncol(S2),ncol(S2)] = (holdS[,,sim-1] %*% holdR[,,sim-1] %*% holdS[,,sim-1])[3,3]

    for(i in 2:(ncol(S2))){
      for(j in 1:(i-1)){
        S2[i,j] = S2[j,i]
      }}

    Sigma = diag(c(0, rep(holdSigE[sim]^2, 2)))

    Sig12 = t(S2[1:(2*r + 1), (2*r + 2):(2*r + 3)])
    Sig22 = S2[1:(2*r + 1), 1:(2*r + 1)]
    Sig11 = S2[(2*r + 2):(2*r + 3), (2*r + 2):(2*r + 3)]

    for(i in 1:n){

      if(trt[i] == 0){

        Sig12 = (S2[ncol(S2) - 1, seq(2, 2*r + 1, 2)])
        Sig11 = S2[ncol(S2) - 1, ncol(S2) - 1]
        Sig22 = S2[seq(2, 2*r + 1, 2), seq(2, 2*r + 1, 2)]

        ymis_mu = (Sig12 %*% ginv(Sig22)) %*% (as.matrix(ST[i, seq(2, 2*r + 1, 2)] - cbind(xbeta_all2)[i, seq(2, 2*r + 1, 2)]))
        ymis_sig = Sig11 - (Sig12) %*% ginv(Sig22) %*% (Sig12)
        ymis = mvrnorm(1, ymis_mu, (ymis_sig))
        ymisb01[i] = ymis

      }

      Zs = rbind(c(1,rep(0, r)), c(0,rep(1, r)))

      if(trt[i] == 1){
        Sig11 = S2[ncol(S2), ncol(S2)]
        Sig22 = S2[1 + seq(0, 2*r + 1, 2),1 + seq(0, 2*r + 1, 2)]
        Sig12 = S2[ncol(S2),1 + seq(0, 2*r + 1, 2)]

        ymis_mu = (Sig12 %*% ginv(Sig22)) %*% (as.matrix(ST[i, 1 + seq(0, 2*r + 1, 2)] - cbind(xbeta_all2)[i,1 + seq(0, 2*r + 1, 2)]))
        ymis_sig = Sig11 - (Sig12) %*% ginv(Sig22) %*% (Sig12)
        ymis = mvrnorm(1, ymis_mu, (ymis_sig))
        ymisb11[i] = ymis
      }
    }

    tmp2 = ymisb01[trt == 0]; tmp3 = ymisb11[trt == 1]

    ST[,2*r + 2] = ymisb01; ST[,2*r + 3] = ymisb11

    x = matrix(c(rep(NA, 25)), ncol = 5, byrow = F)
    summand = rep(list(x), n/2)
    x = matrix(c(rep(NA, 5)), ncol = 1, byrow = F)
    summand2 = rep(list(x), n/2)

    # estimate mean parameters
    x = matrix(c(rep(1, r)), ncol = 1, byrow = F)

    X1 = rep(list(x), n)

    for(i in 1:(n/2)){
      summand[[i]] = t((X1[[i]])) %*% (X1[[i]])
      summand2[[i]] = t(X1[[i]]) %*% ST[i,seq(2, 2*r + 1, 2)] - c(ST[i,ncol(ST) - 1]) # # update to make more general
    }

    betat0 = Reduce(`+`, lapply(summand2, function(x) replace(x, is.na(x), 0)))/n*2/3

    V = X1
    for(i in 1:n){
      V[[i]] = ((X1[[i]]) %*% t(X1[[i]]))[1,1]
    }

    V = ginv(Reduce("+", V))
    betas = mvrnorm(1,(betat0), holdSigE[sim] * V)
    holdbeta0[sim] = betas[1];

    summand2 = rep(list(x), n/2)

    for(i in (n/2+1):n){
      summand[[i - n/2]] = t((X1[[i]])) %*% (X1[[i]])
      summand2[[i - (n/2)]] = t(X1[[i]]) %*% ST[i,1 + seq(2, 2*r + 1, 2)] - c(ST[i,ncol(ST)])
    }

    betat1 = Reduce(`+`, lapply(summand2, function(x) replace(x, is.na(x), 0))) / n * 2 / 3

    V = X1
    for(i in 1:n){
      V[[i]] = ((X1[[i]]) %*% t(X1[[i]]))[1,1]
    }

    V = ginv(Reduce("+", V))
    betas = mvrnorm(1, (betat1), holdSigE[sim] * V)

    holdbeta1[sim] = betas[1]

    Xmat = rep(1, n/2); obsS1 = ST[trt==1,1]
    Lambda0t = matrix(c(rep(0.1, 1)), nrow = 1); tauST0 = holdS[1,1,sim-1]^2
    v = ginv(Lambda0t + as.numeric(tauST0) * (t(Xmat) %*% Xmat))
    m = v %*% (tauST0*t(Xmat) %*% as.matrix(obsS1))
    betaT=c(rmvnorm(1,m,v/n))

    holdalpha1[sim] = betaT[1]

    tmp1 = obsS1[1:(n/2)] - Xmat[1:(n/2)] * betaT

    holdmu[,sim] = c(holdalpha1[sim], holdbeta0[sim], holdbeta1[sim])

    # estimate sd parameters
    s1 =  rinvgamma(1, shape = a + n/4, scale = (sum(tmp1^2)/2 + b))
    s2 =  rinvgamma(1, shape = a + n/4, scale = (sum(tmp2^2)/2 + b))
    s3 =  rinvgamma(1, shape = a + n/4, scale = (sum(tmp3^2)/2 + b))
    holdS[,,sim] = S = diag(c(sqrt(s1), sqrt(s2), sqrt(s3)))

    r23 = holdR[2,3,sim-1]; r13 = holdR[1,3,sim]; r12 = holdR[1,2,sim-1]
    resid = cbind(tmp1, tmp3)

    # estimate correlation parameters
    phat = cor(tmp1, tmp3)

    y = rnorm(1, 0.5 * log((1+holdR[1,3,sim-1]) / (1- holdR[1,3,sim-1])), sd = sqrt(1/(n-3)))

    R2 = holdR[,,sim-1]; R2[1,3] = R2[3,1] = ifisherz(y)
    if(condindfit) R2[1,2] = R2[2,1] = R2[1,3] * R2[2,3]

    summand = apply(resid, 1, function(resid) t(resid) %*% ginv(S[c(1,3), c(1,3)] %*% R2[c(1,3), c(1,3)] %*% S[c(1,3), c(1,3)]) %*% (resid) )

    ratio = exp(fdelt(n/2, R = R2, j = sum(summand)))*(1/(1- holdR[1,3,sim-1]^2))

    R2 = holdR[,,sim-1]; R2[1,3] = R2[3,1] = holdR[1,3,sim-1]
    if(condindfit){ R2[1,2] = R2[2,1] = R2[1,3] * R2[2,3] }

    summand = apply(resid, 1, function(resid) t(resid) %*% ginv(S[c(1,3), c(1,3)] %*% R2[c(1,3), c(1,3)] %*% S[c(1,3), c(1,3)]) %*% (resid) )

    ratio2 = exp(fdelt(n/2, R = R2, j = sum(summand)))*(1/(1-ifisherz(y)^2))

    prob = max(0, min(1, (ratio / ratio2)))
    if(is.na(prob)) next

    z = rbern(1, prob)
    if(z == 1 & sim > burnin) accept = accept + 1

    R = holdR[,,sim-1]
    r13 = ifisherz(y)*z + (1-z) * holdR[1,3,sim-1]
    r13 = r13[1]
    R[1,3] = R[3,1] = r13

    R[2,3] = R[3,2] = holdR[2,3,sim] = rBeta_ab(1, shape1 = 8, shape2 = 5, a = -1, b = 1)


    if(prior == 'unif'){R[2,3] = R[3,2] = holdR[2,3,sim] = runif(1,-1,1)}

    R[1,2] = R[2,1] = runif(1,-1,1)
    if(condindfit){ R[1,2] = R[2,1] = R[1,3] * R[2,3] }


    if(any(eigen(R)$values<0)) next; if(any(abs(R)>1)) next

    holdR[,, sim] = R; holdS[,,sim] = S

    # calculate slope and intercept terms: gamma_1 and gamma_0
    slope[sim] = (holdR[1,3,sim] * (holdS[3,3,sim]) - holdR[1,2,sim] * (holdS[2,2,sim])) / (holdS[1,1,sim])
    int[sim] = (holdmu[3,sim]- holdmu[2,sim]) - slope[sim]*(holdmu[1,sim])

    sim = sim + 1
    if(sim %% 50 == 0) print(sim)

  }


  params = data.frame(int = int, slope = slope, r13 = holdR[1,3,], r12 = holdR[1,2,], r23 = holdR[2,3,],
                      s1 = holdS[1,1,], s2 = holdS[2,2,], s3 = holdS[3,3,],
                      beta0 = holdbeta0, beta1 = holdbeta1, alpha1 = holdalpha1,
                      omega1 = holdomega1, omega2 = holdomega2, psi2 = holdpsi2,
                      psi1 = holdpsi2, b0 = holdB0, b1 = holdB1
  )

  result = list(params = params,
                accept.rate = accept/(SIM - burnin),
                args = list(SIM = SIM, burnin = burnin, n = n))


  return(result)

}
