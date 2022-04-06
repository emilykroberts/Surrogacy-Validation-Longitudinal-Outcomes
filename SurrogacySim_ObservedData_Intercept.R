# load required packages
library(corpcor); library(bayesSurv); library(MASS); library(GenKern); library(coda); library(psych)
library(mvtnorm); library(rootSolve); library(Matrix); library(MCMCpack); library(LaplacesDemon)
library(reshape2); library(tidyr); library(ggplot2); library(pracma); library(lme4); library(lcmm)
library(HardyWeinberg); library(ExtDist); library(lcmix)

# set up needed simulation settings to run on cluster -> comment out as needed
set.seed(323);
array_id2 = array_id = 1; array_id2 = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(323 + array_id2);
sim = 2

# simulation parameters
burnin = 1; n = 300; SIM = 10000 # burnin, sample size, and number of MCMC iterations
accept = 0 # counter for acceptance ratio
condindfit = T # fit conditional independence, boolean
prior = 'beta' # or unif
misspecify = T

# generating values
alpha1 = 2; psi2 = 0
beta0 = 2; omega1 = 0
beta1 = 3.1; omega2 = 0
eps1 = 1; eps2 = 1; eps3 = 1; tau4 = 0.5; eta4 = 1
theta10 = 0.15; theta11 = 0.7; thetaT = 0.2142857
sigmae = 0.3; tau = sigmae^2 # represents ei
R = matrix(rep(1,3*3),3,3); R[1,2] = R[2,1] = theta10; R[1,3] = R[3,1] =theta11; R[2,3] = R[3,2] = thetaT
S = diag(c(0.25, 0.25, 0.25)); 

### generate longitudinal T0 and T1 = mu + ui + eij #ui shared across outcomes per person
trueR = biR = R
mu = c(alpha1, beta0, beta1) 

Xi = cbind(rep(1, n))

allsamp = mvrnorm(n, mu = c(alpha1, beta0, beta1, beta0, beta1, 
                            beta0, beta1) , diag(c(tau, tau, tau, tau, tau, tau, tau)))

allsamp[,1] = mu[1]

# allow for misspecifying MVN distribution of random effects
uitruesamp = uitrue_normal = uitrue = mvrnorm(n, c(0, 0, 0), S %*% biR %*% S)

if(misspecify){
  uitrue = rmvgamma(n, shape = 0.5, rate = 0.5, corr = S %*% biR %*% S)
  uitruesamp = uitrue = uitrue - colMeans(uitrue)
  # hist(uitrue[,1]); hist(uitrue[,2]); hist(uitrue[,3]) 
}

allsamp[, 1:3] = allsamp[, 1:3] + uitrue
allsamp[, 4:5] = allsamp[, 4:5] + uitrue[,2:3]
allsamp[, 6:7] = allsamp[, 6:7] + uitrue[,2:3]

# needed hyperparameters
a = b = 0.1

# set up matrices to save results
# save data corresponding to this simulation on cluster

trt = c(rep(0,n)); trt[(1*n-(n/2-1)):(n*1)] = 1
ST = samp = allsamp

holdmu = matrix(rep(0,3*SIM),3,SIM); holdS = array(rep(0,3*3*SIM),dim=c(3,3,SIM))
holdR = array(rep(0,3*3*SIM),dim = c(3,3,SIM)); holdR[,,1] = R; holdS[,,1] = S
holdpsi1 = (rep(0,1*SIM)); holdpsi2 = (rep(0,1*SIM))
holdomega1 = (rep(0,1*SIM)); holdomega2 = (rep(0,1*SIM))
holdalpha0 = (rep(0,1*SIM)); holdalpha1 = (rep(0,1*SIM)); 
holdbeta0 = (rep(0,1*SIM)); holdbeta1 = (rep(0,1*SIM)); 

slope = array(0,c((sim-2),1))
intmarginalize = (rep(0,1*SIM)); 
slopemarginalize = (rep(0,1*SIM)); 
intmarginalizeclosed = (rep(0,1*SIM)); 
slopemarginalizeclosed = (rep(0,1*SIM)); 

holdB0 = array(0,c((sim-2),1))
holdB1 = array(0,c((sim-2),1))
holdSigE = array(0,c((sim-2),1))
holdB0[1] = holdB1[1] = sqrt(sigmae)
int = array(0,c((sim-1),1))
int1 = array(0,c((sim-1),1))

holdalpha1[1] = mu[1]; holdbeta0[1] = mu[2]; holdbeta1[1] = mu[3]
bivec = uitruesamp[,2]*(1-trt) + uitruesamp[,3]*trt

r=3

x_1 = z_1 = matrix(c(c(1, 1, 1, 1)), ncol = 2, byrow = F)
Zi = Zi0 = Zi1 = rep(list(z_1), n)
X1 = rep(list(x_1, n))

for(i in 1:n){
  
  X1[[i]] = matrix(c(c(1, 1, 1, 1, 1, 1)), ncol = 2, byrow = T)
  Zi1[[i]] = matrix(c(c(1, 1, 1)), ncol = 3, byrow = T)
  Zi[[i]] = matrix(c(c(1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1)), ncol = 2, byrow = T)
  
}

ST = cbind(ST, uitruesamp[,2:3])

# likelihood calculation for MH step
fdelt = function(n,R,j){ return(-(n/2)*log(det(R))+(-0.5*j) )}

# start simulation

while(sim <= SIM){
  
  mu = cbind(holdalpha1[sim-1], holdbeta0[sim-1], holdbeta1[sim-1], holdbeta0[sim-1], holdbeta1[sim-1])
  
  # estimate error term
  Sig = holdS[,,sim-1] %*% holdR[,,sim-1] %*% holdS[,,sim-1]
  tmp5 = c(ST[trt == 1,3] - Xi[trt == 1] * holdbeta1[sim-1] - 1 * c(bivec[trt == 1]))
  bivec[trt == 1] == uitrue[trt == 1,3]
  tmp4 = c(ST[trt == 0,2] - Xi[trt == 0] * holdbeta0[sim-1] - 1 * c(bivec[trt == 0]))
  
  B0 = B1 = rinvgamma(1, shape = a + length(c(tmp4, tmp5))/2, scale = (sum((tmp4 + tmp5)^2)/2 + b)) 
  
  holdB0[sim] = B0; holdB1[sim] = B1
  holdSigE[sim] = sqrt(B0)
  Sigma = diag(c(rep(holdSigE[sim]^2, 3)))
  
  # estimate random effects
  xbeta_all2 = ST
  xbeta_all2[,1] = holdalpha1[sim-1]
  xbeta_all2[,c(2,4,6)] = holdbeta0[sim-1]
  xbeta_all2[,c(3,5,7)] = holdbeta1[sim-1]
  ymisb01 = ymisb11 = ymiss = rep(NA, n)
  
  S2 = cbind(holdS[,,sim-1] %*% holdR[,,sim-1] %*% holdS[,,sim-1], 0, 0, 0, 0, 0, 0)
  S2 = rbind(S2, 0, 0, 0, 0, 0, 0)
  S2[1, 2:4] = (holdS[,,sim-1] %*% holdR[,,sim-1] %*% holdS[,,sim-1])[1,2] # cov(S1, T0)
  S2[1, 5:7] = (holdS[,,sim-1] %*% holdR[,,sim-1] %*% holdS[,,sim-1])[1,3] # cov(S1, T1)
  S2[2:4, 2:4] = (holdS[,,sim-1] %*% holdR[,,sim-1] %*% holdS[,,sim-1])[2,2] 
  S2[5:7, 5:7] = (holdS[,,sim-1] %*% holdR[,,sim-1] %*% holdS[,,sim-1])[3,3] 
  S2[2:4, 5:7] = (holdS[,,sim-1] %*% holdR[,,sim-1] %*% holdS[,,sim-1])[2,3] # cov(T0, T1) = cov(b0, b1)
  S2[3,3] = S2[4,4] = S2[2,2] = (holdS[,,sim-1] %*% holdR[,,sim-1] %*% holdS[,,sim-1])[2,2] + holdSigE[sim]^2
  S2[5,5] = S2[6,6] = S2[7,7] = (holdS[,,sim-1] %*% holdR[,,sim-1] %*% holdS[,,sim-1])[3,3] + holdSigE[sim]^2
  S2[1,8] = (holdS[,,sim-1] %*% holdR[,,sim-1] %*% holdS[,,sim-1])[1,2]
  S2[1,9] = (holdS[,,sim-1] %*% holdR[,,sim-1] %*% holdS[,,sim-1])[1,3]
  S2[2,8] = S2[3,8] = S2[4,8] = (holdS[,,sim-1] %*% holdR[,,sim-1] %*% holdS[,,sim-1])[2,2]
  S2[2,9] = S2[3,9] = S2[4,9] = (holdS[,,sim-1] %*% holdR[,,sim-1] %*% holdS[,,sim-1])[2,3]
  
  S2[5,8] = S2[6,8] = S2[7,8] = (holdS[,,sim-1] %*% holdR[,,sim-1] %*% holdS[,,sim-1])[2,3]
  S2[5,9] = S2[6,9] = S2[7,9] = (holdS[,,sim-1] %*% holdR[,,sim-1] %*% holdS[,,sim-1])[3,3]
  S2[8,9] = (holdS[,,sim-1] %*% holdR[,,sim-1] %*% holdS[,,sim-1])[2,3] # b01 b02
  S2[8,8] = (holdS[,,sim-1] %*% holdR[,,sim-1] %*% holdS[,,sim-1])[2,2]
  S2[9,9] = (holdS[,,sim-1] %*% holdR[,,sim-1] %*% holdS[,,sim-1])[3,3]
  
  for(i in 2:(9)){
    for(j in 1:(i-1)){
      S2[i,j] = S2[j,i]
    }}
  
  Sigma = diag(c(0, rep(holdSigE[sim]^2, 2)))
  
  Sig12 = t(S2[c(1:7), c(8:9)])
  Sig22 = S2[c(1:7), c(1:7)]
  Sig11 = S2[c(8:9), c(8:9)]
  
  for(i in 1:n){
    
    if(trt[i] == 0){
      
      Sig12 = (S2[c(8), c(2, 4, 6)])
      Sig11 = S2[c(8), c(8)]
      Sig22 = S2[c(2,4,6), c(2,4,6)]
      
      ymis_mu = (Sig12 %*% ginv(Sig22)) %*% (as.matrix(ST[i, c(2,4,6)] - cbind(xbeta_all2)[i, c(2,4,6)])) 
      ymis_sig = Sig11 - (Sig12) %*% ginv(Sig22) %*% (Sig12) 
      ymis = mvrnorm(1, ymis_mu, (ymis_sig))
      ymisb01[i] = ymis
      
    }
    
    Zs = rbind(c(1,0,0,0), c(0,1,1,1))
    
    if(trt[i] == 1){
      Sig11 = S2[c(9), c(9)]
      Sig22 = S2[c(1, 3, 5, 7),c(1, 3, 5, 7)]
      Sig12 = S2[c(9),c(1, 3, 5, 7)]
      
      ymis_mu = (Sig12 %*% ginv(Sig22)) %*% (as.matrix(ST[i, c(1, 3, 5, 7)] - cbind(xbeta_all2)[i,c(1, 3, 5, 7)])) 
      ymis_sig = Sig11 - (Sig12) %*% ginv(Sig22) %*% (Sig12) 
      ymis = mvrnorm(1, ymis_mu, (ymis_sig))
      ymisb11[i] = ymis
    }
  }
  
  tmp1 = ST[trt == 1,1]; tmp2 = ymisb01[trt == 0]; tmp3 = ymisb11[trt == 1]
  
  ST[,8] = ymisb01; ST[,9] = ymisb11
  
  x = matrix(c(rep(NA, 25)), ncol = 5, byrow = F)
  summand = rep(list(x), n/2)
  x = matrix(c(rep(NA, 5)), ncol = 1, byrow = F)
  summand2 = rep(list(x), n/2)
  
  r = 3
  
  # estimate mean parameters
  x = matrix(c(rep(1, r)), ncol = 1, byrow = F)
  
  X1 = rep(list(x), n)
  
  for(i in 1:(n/2)){
    summand[[i]] = t((X1[[i]])) %*% (X1[[i]])
    summand2[[i]] = t(X1[[i]]) %*% ST[i,c(2,4,6)] - c(ST[i,c(8)])
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
    summand2[[i - (n/2)]] = t(X1[[i]]) %*% ST[i,c(3,5,7)] - c(ST[i,c(9)])
  }
  
  betat1 = Reduce(`+`, lapply(summand2, function(x) replace(x, is.na(x), 0))) / n * 2 / 3
  
  V = X1
  for(i in 1:n){ 
    V[[i]] = ((X1[[i]]) %*% t(X1[[i]]))[1,1]
  }
  
  V = ginv(Reduce("+", V))
  betas = mvrnorm(1, (betat1), holdSigE[sim] * V) 
  
  holdbeta1[sim] = betas[1]
  
  Xmat = Xi[trt==1]; obsS1 = ST[trt==1,1]
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
  holdS[,,sim] = Sigma = S = diag(c(sqrt(s1), sqrt(s2), sqrt(s3)))
  
  r23 = holdR[2,3,sim-1]; r13 = holdR[1,3,sim]; r12 = holdR[1,2,sim-1]
  resid = cbind(tmp1, tmp3)
  
  # estimate correlation parameters
  if(T){
    phat = cor(tmp1, tmp3)
    
    y = rnorm(1, 0.5 * log((1+holdR[1,3,sim-1]) / (1- holdR[1,3,sim-1])), sd = sqrt(1/(n-3)))
    
    R2 = holdR[,,sim-1]; R2[1,3] = R2[3,1] = ifisherz(y)
    if(condindfit) R2[1,2] = R2[2,1] = R2[1,3] * R2[2,3]
    
    resid1 = (0.5 * log((1 + holdR[1,3,sim-1]) / (1 - holdR[1,3,sim-1]))) - y
    
    summand = apply(resid, 1, function(resid) t(resid) %*% ginv(S[c(1,3), c(1,3)] %*% R2[c(1,3), c(1,3)] %*% S[c(1,3), c(1,3)]) %*% (resid) )
    summand1 = resid1^2/(n-3)
    
    ratio = exp(fdelt(n/2, R = R2, j = sum(summand)))*(1/(1- holdR[1,3,sim-1]^2))
    
    R2 = holdR[,,sim-1]; R2[1,3] = R2[3,1] = holdR[1,3,sim-1]
    if(condindfit){ R2[1,2] = R2[2,1] = R2[1,3] * R2[2,3] }
    
    resid1 = y - (0.5 * log((1 + holdR[1,3,sim-1]) / (1- holdR[1,3,sim-1])))
    
    summand = apply(resid, 1, function(resid) t(resid) %*% ginv(S[c(1,3), c(1,3)] %*% R2[c(1,3), c(1,3)] %*% S[c(1,3), c(1,3)]) %*% (resid) )
    summand1 = resid1^2 / (n-3)
    
    ratio2 = exp(fdelt(n/2, R = R2, j = sum(summand)))*(1/(1-ifisherz(y)^2))
    
    prob = max(0, min(1, (ratio / ratio2)))
    if(is.na(prob)) next
    
    z = rbern(1, prob); if(z == 1) accept = accept + 1
    
    r13 = ifisherz(y)*z + (1-z) * holdR[1,3,sim-1]
    r13 = r13[1]
    R[1,3] = R[3,1] = r13
    
    R[2,3] = R[3,2] = holdR[2,3,sim] = rBeta_ab(1, shape1 = 8, shape2 = 5, a = -1, b = 1)
    
    if(prior == 'unif'){R[2,3] = R[3,2] = holdR[2,3,sim] = runif(1,-1,1)}
    
    R[1,2] = R[2,1] = runif(1,-1,1)
    if(condindfit){ R[1,2] = R[2,1] = R[1,3] * R[2,3] }
  }
  
  if(any(eigen(R)$values<0)) next; if(any(abs(R)>1)) next
  
  holdR[,, sim] = R; holdS[,,sim] = S
  
  # calculate slope and intercept terms: gamma_1 and gamma_0
  slope[sim] = (holdR[1,3,sim] * (holdS[3,3,sim]) - holdR[1,2,sim] * (holdS[2,2,sim])) / (holdS[1,1,sim])
  
  int[sim] = (holdmu[3,sim]- holdmu[2,sim]) - slope[sim]*(holdmu[1,sim])
  
  # print(R)
  # print(S)
  
  sim = sim + 1
  print(sim)
  
}  

### traceplots
plot(holdR[1,2,])
plot(holdR[2,3,])
plot(holdR[1,3,])

plot(int)
plot(slope)

## save results
params = data.frame(beta0SE = numeric(1), beta1SE = numeric(1), alpha0SE = numeric(1), alpha1SE = numeric(1),
                    psi1SE = numeric(1), omega1SE = numeric(1), psi2SE = numeric(1), omega2SE = numeric(1),
                    sigs0 = numeric(1), SEsigs0 = numeric(1), sigs1 = numeric(1), SEsigs1 = numeric(1),
                    sigmat0 = numeric(1), SEsigmat0 = numeric(1),sigmat1 = numeric(1), SEsigmat1 = numeric(1), 
                    ps = numeric(1), SEps = numeric(1),p00 = numeric(1), SEp00 = numeric(1),
                    p01 = numeric(1),SEp01 = numeric(1),p10 = numeric(1),SEp10 = numeric(1),p11 = numeric(1),
                    SEp11 = numeric(1),pt = numeric(1),SEpt = numeric(1), S0l = numeric(1), S0u = numeric(1),
                    S1l = numeric(1), S1u = numeric(1), T0l = numeric(1),T0u = numeric(1), T1l = numeric(1),
                    T1u = numeric(1), psl = numeric(1), PSu = numeric(1),p00l = numeric(1), p00u = numeric(1),
                    p01l = numeric(1),p01u = numeric(1),p10l = numeric(1),p10u = numeric(1),
                    r1 = numeric(1), SEr1 = numeric(1), r2= numeric(1), SEr2= numeric(1),r3 = numeric(1),SEr3 = numeric(1),
                    r1l = numeric(1), r1u = numeric(1),r2l = numeric(1), r2u = numeric(1),
                    r3l = numeric(1),r3u = numeric(1), b0 = numeric(1), b0SE = numeric(1), b1 = numeric(1), b1SE = numeric(1))

prentice = data.frame(bet_s = numeric(1), SEbs = numeric(1),bet_trtt = numeric(1), SEbtt = numeric(1), 
                      bet_trts = numeric(1), SEbts = numeric(1), gam_trt = numeric(1), SEgt = numeric(1),
                      gam_s = numeric(1), SEgs = numeric(1))

pren.post = data.frame(bet_sp= numeric(1), SEbet_sp= numeric(1),bet_trttp= numeric(1), SEbet_trttp= numeric(1), 
                       bet_trtsp= numeric(1), SEbet_trtsp= numeric(1), gam_trtp= numeric(1), SEgam_trtp= numeric(1), 
                       gam_sp= numeric(1),SEgam_sp= numeric(1))

PS = data.frame(dat_int = numeric(1), dat_intSE = numeric(1),dat_sl = numeric(1), dat_slSE = numeric(1), 
                mean_int = numeric(1), SEmean_int = numeric(1),L_int = numeric(1), U_int = numeric(1),mean_sl = numeric(1)
                ,SEmean_sl = numeric(1),L_sl = numeric(1),U_sl = numeric(1),
                postdat_int = numeric(1), postdat_intSE = numeric(1), postdat_sl = numeric(1), postdat_slSE = numeric(1),
                mean_int = numeric(1), SEmean_int = numeric(1),L_int = numeric(1), U_int = numeric(1), covslint = numeric(1),covslint1 = numeric(1),
                int_coveragE = numeric(1), slope_coveragE = numeric(1), int1_coveragE = numeric(1),
                intmarg_mean = numeric(1), intmarg_SE = numeric(1), intmargl = numeric(1), intmargu = numeric(1), intmarg_coveragE = numeric(1),
                intaccept = numeric(1), int1accept = numeric(1), intmargaccept = numeric(1), slopeaccept = numeric(1),
                slopemarg_mean = numeric(1), slopemarg_SE = numeric(1), slopemargl = numeric(1), slopemargu = numeric(1),
                intmarg_meancl = numeric(1), intmarg_SEcl = numeric(1), intmarglcl = numeric(1), intmargucl = numeric(1),
                slopemarg_meancd = numeric(1), slopemarg_SEcd = numeric(1), slopemarglcd = numeric(1), slopemargucd = numeric(1))

covs = data.frame(psl = numeric(1), psu = numeric(1),p00l = numeric(1), p00u = numeric(1),
                  p01l = numeric(1),p01u = numeric(1),p10l = numeric(1),p10u = numeric(1),p11l = numeric(1),
                  p11u = numeric(1),ptl = numeric(1),ptu = numeric(1),s0l = numeric(1),s0u = numeric(1)
                  ,s1l = numeric(1),s1u = numeric(1),t0l = numeric(1),t0u = numeric(1),t1l = numeric(1),t1u = numeric(1),
                  psind = numeric(1),p00ind = numeric(1),p01ind = numeric(1),p10ind = numeric(1),p11ind = numeric(1),ptind = numeric(1),
                  s0ind = numeric(1),s1ind = numeric(1),t0ind = numeric(1),t1ind = numeric(1))

params[1] = sqrt(var(holdbeta0[burnin:(sim-1)]))
params[2] = sqrt(var(holdbeta1[burnin:(sim-1)]))
#params[3] = sqrt(var(holdalpha0[burnin:(sim-1)]))
params[4] = sqrt(var(holdalpha1[burnin:(sim-1)]))
params[5] = sqrt(var(holdpsi1[burnin:(sim-1)]))
params[6] = sqrt(var(holdomega1[burnin:(sim-1)]))
params[7] = sqrt(var(holdpsi2[burnin:(sim-1)]))
params[8] = sqrt(var(holdomega2[burnin:(sim-1)]))

params[9] = mean(holdS[1,1,burnin:(sim-1)])
params[10] = sqrt(var(holdS[1,1,burnin:(sim-1)]))
params[11] = mean(holdS[2,2,burnin:(sim-1)])
params[12] = sqrt(var(holdS[2,2,burnin:(sim-1)]))
params[13] = mean(holdS[3,3,burnin:(sim-1)])
params[14] = sqrt(var(holdS[3,3,burnin:(sim-1)]))

params[17] = mean(holdR[1,2,burnin:(sim-1)])
params[18] = sqrt(var(holdR[1,2,burnin:(sim-1)]))
params[19] = mean(holdR[1,3,burnin:(sim-1)])
params[20] = sqrt(var(holdR[1,3,burnin:(sim-1)]))
params[23] = mean(holdR[2,3,burnin:(sim-1)])
params[24] = sqrt(var(holdR[2,3,burnin:(sim-1)]))

params[29] = quantile(holdS[1,1,burnin:(sim-1)], probs = 0.025, na.rm = TRUE)
params[30] = quantile(holdS[1,1,burnin:(sim-1)], probs = 0.975, na.rm = TRUE)
params[31] = quantile(holdS[2,2,burnin:(sim-1)], probs = 0.025, na.rm = TRUE)
params[32] = quantile(holdS[2,2,burnin:(sim-1)], probs = 0.975, na.rm = TRUE)
params[33] = quantile(holdS[3,3,burnin:(sim-1)], probs = 0.025, na.rm = TRUE)
params[34] = quantile(holdS[3,3,burnin:(sim-1)], probs = 0.975, na.rm = TRUE)

params[37] = quantile(holdR[1,2,burnin:(sim-1)], probs = 0.025, na.rm = TRUE)
params[38] = quantile(holdR[1,2,burnin:(sim-1)], probs = 0.975, na.rm = TRUE)
params[39] = quantile(holdR[1,3,burnin:(sim-1)], probs = 0.025, na.rm = TRUE)
params[40] = quantile(holdR[1,3,burnin:(sim-1)], probs = 0.975, na.rm = TRUE)
params[43] = quantile(holdR[2,3,burnin:(sim-1)], probs = 0.025, na.rm = TRUE)
params[44] = quantile(holdR[2,3,burnin:(sim-1)], probs = 0.975, na.rm = TRUE)

params[57] = mean(holdB0[burnin:(sim-1)])
params[58] = sqrt(var(holdB0[burnin:(sim-1)]))
params[59] = mean(holdB1[burnin:(sim-1)])
params[60] = sqrt(var(holdB1[burnin:(sim-1)]))


covs[1] = quantile(holdR[1,2,burnin:sim-1], probs = 0.025, na.rm = T)
covs[2] = quantile(holdR[1,2,burnin:sim-1], probs = 0.975, na.rm = T)
covs[3] = quantile(holdR[1,3,burnin:sim-1], probs = 0.025, na.rm = T)
covs[4] = quantile(holdR[1,3,burnin:sim-1], probs = 0.975, na.rm = T)
covs[7] = quantile(holdR[2,3,burnin:sim-1], probs = 0.025, na.rm = T)
covs[8] = quantile(holdR[2,3,burnin:sim-1], probs = 0.975, na.rm = T)

covs[13] = quantile(holdS[1,1,burnin:sim-1], probs = 0.025, na.rm = T)
covs[14] = quantile(holdS[1,1,burnin:sim-1], probs = 0.975, na.rm = T)
covs[15] = quantile(holdS[2,2,burnin:sim-1], probs = 0.025, na.rm = T)
covs[16] = quantile(holdS[2,2,burnin:sim-1], probs = 0.975, na.rm = T)
covs[17] = quantile(holdS[3,3,burnin:sim-1], probs = 0.025, na.rm = T)
covs[18] = quantile(holdS[3,3,burnin:sim-1], probs = 0.975, na.rm = T)

covs[21] = as.numeric((holdR[1,2,1]>covs[1])&(holdR[1,2,1]<covs[2]))
covs[22] = as.numeric((holdR[1,3,1]>covs[3])&(holdR[1,2,1]<covs[4]))
covs[24] = as.numeric((holdR[2,3,1]>covs[7])&(holdR[2,3,1]<covs[8]))
covs[27] = as.numeric((holdS[1,1,1]>covs[13])&(holdS[1,1,1]<covs[14]))
covs[28] = as.numeric((holdS[2,2,1]>covs[15])&(holdS[2,2,1]<covs[16]))

PS[5] = mean(int[burnin:sim-1], na.rm = T)
PS[6] = sqrt(var(int[burnin:sim-1], na.rm = T))
PS[7] = quantile(int[burnin:sim-1], probs = 0.025, na.rm = T)
PS[8] = quantile(int[burnin:sim-1], probs = 0.975, na.rm = T)
PS[9] = mean(slope[burnin:sim-1], na.rm = T)
PS[10] = sqrt(var(slope[burnin:sim-1], na.rm = T))
PS[11] = quantile(slope[burnin:sim-1], probs = 0.025, na.rm = T)
PS[12] = quantile(slope[burnin:sim-1], probs = 0.975, na.rm = T)

PS[17] = mean(int1[burnin:sim-1], na.rm = T)
PS[18] = sqrt(var(int1[burnin:sim-1], na.rm = T))
PS[19] = quantile(int1[burnin:sim-1], probs = 0.025, na.rm = T)
PS[20] = quantile(int1[burnin:sim-1], probs = 0.975, na.rm = T)
PS[21] = cov(int[burnin:sim-1],slope[burnin:sim-1])
PS[22] = cov(int1[burnin:sim-1],slope[burnin:sim-1])
PS[22] = cov(int1[burnin:sim-1],slope[burnin:sim-1])
PS[23] = as.numeric(PS[7]<int[1]&PS[8]>int[1])
PS[24] = as.numeric(PS[11]<slope[1]&PS[12]>slope[1])
PS[25] = as.numeric(PS[19]<int1[1]&PS[20]>int1[1])
PS[26] = mean(intmarginalize[burnin:sim-1], na.rm = T)
PS[27] = sqrt(var(intmarginalize[burnin:sim-1], na.rm = T))
PS[28] = quantile(intmarginalize[burnin:sim-1], probs = 0.025, na.rm = T)
PS[29] = quantile(intmarginalize[burnin:sim-1], probs = 0.975, na.rm = T)
PS[30] = as.numeric(PS[28]< 0 & PS[29]> 0) # may have to fix
PS[31] = as.numeric(PS[7]<0 & PS[8]>=0)
PS[32] = as.numeric(PS[19]<0 & PS[20]>=0)
PS[33] = as.numeric(PS[28]<0 & PS[29]>=0)
PS[34] = as.numeric(PS[11]>0)
PS[35] = mean(slopemarginalize[burnin:sim-1], na.rm = T)
PS[36] = sqrt(var(slopemarginalize[burnin:sim-1], na.rm = T))
PS[37] = quantile(slopemarginalize[burnin:sim-1], probs = 0.025, na.rm = T)
PS[38] = quantile(slopemarginalize[burnin:sim-1], probs = 0.975, na.rm = T)
PS[39] = mean(intmarginalizeclosed[burnin:sim-1], na.rm = T)
PS[40] = sqrt(var(intmarginalizeclosed[burnin:sim-1], na.rm = T))
PS[41] = quantile(intmarginalizeclosed[burnin:sim-1], probs = 0.025, na.rm = T)
PS[42] = quantile(intmarginalizeclosed[burnin:sim-1], probs = 0.975, na.rm = T)
PS[43] = mean(slopemarginalizeclosed[burnin:sim-1], na.rm = T)
PS[44] = sqrt(var(slopemarginalizeclosed[burnin:sim-1], na.rm = T))
PS[45] = quantile(slopemarginalizeclosed[burnin:sim-1], probs = 0.025, na.rm = T)
PS[46] = quantile(slopemarginalizeclosed[burnin:sim-1], probs = 0.975, na.rm = T)

estcoef = data.frame(beta0mean = numeric(1), beta1mean = numeric(1), omega1mean = numeric(1), alpha0mean = numeric(1),
                     alpha1mean = numeric(1),psi1mean = numeric(1),psi2mean = numeric(1),omega2mean = numeric(1),
                     beta0l = numeric(1),beta0u = numeric(1), beta1l = numeric(1),beta1u = numeric(1),omega1l = numeric(1),
                     omega1u = numeric(1), alpha0l = numeric(1),alpha0u = numeric(1),
                     alpha1l = numeric(1),alpha1u = numeric(1),psi1l = numeric(1),psi1u = numeric(1),
                     psi2l = numeric(1),psi2u = numeric(1),omega2l = numeric(1),omega2u = numeric(1), 
                     beta0ind = numeric(1), beta1ind = numeric(1),omega1ind = numeric(1), alpha0ind = numeric(1),
                     alpha1ind = numeric(1), psi1ind = numeric(1),psi2ind = numeric(1),omega2ind = numeric(1))

estcoef[1] = mean(holdbeta0[burnin:sim-1], na.rm = T)
estcoef[2] = mean(holdbeta1[burnin:sim-1], na.rm = T)
estcoef[3] = mean(holdomega1[burnin:sim-1], na.rm = T)
#estcoef[4] = alpha0mean = mean(holdalpha0[burnin:sim-1], na.rm = T)
estcoef[5] = mean(holdalpha1[burnin:sim-1], na.rm = T)
estcoef[6] = mean(holdpsi1[burnin:sim-1], na.rm = T)
estcoef[7] = mean(holdpsi2[burnin:sim-1], na.rm = T)
estcoef[8] = omega2mean = mean(holdomega2[burnin:sim-1], na.rm = T)
estcoef[9] = quantile(holdbeta0[burnin:sim-1], probs = 0.025, na.rm = T)
estcoef[10] = quantile(holdbeta0[burnin:sim-1], probs = 0.975, na.rm = T)
estcoef[11] = quantile(holdbeta1[burnin:sim-1], probs = 0.025, na.rm = T)
estcoef[12] = quantile(holdbeta1[burnin:sim-1], probs = 0.975, na.rm = T)
estcoef[13] = quantile(holdomega1[burnin:sim-1], probs = 0.025, na.rm = T)
estcoef[14] = quantile(holdomega1[burnin:sim-1], probs = 0.975, na.rm = T)
#estcoef[15] = quantile(holdalpha0[burnin:sim-1], probs = 0.025, na.rm = T)
#estcoef[16] = quantile(holdalpha0[burnin:sim-1], probs = 0.975, na.rm = T)
estcoef[17] = quantile(holdalpha1[burnin:sim-1], probs = 0.025, na.rm = T)
estcoef[18] = quantile(holdalpha1[burnin:sim-1], probs = 0.975, na.rm = T)
estcoef[19] = quantile(holdpsi1[burnin:sim-1], probs = 0.025, na.rm = T)
estcoef[20] = quantile(holdpsi1[burnin:sim-1], probs = 0.975, na.rm = T)
estcoef[21] = quantile(holdpsi2[burnin:sim-1], probs = 0.025, na.rm = T)
estcoef[22] = quantile(holdpsi2[burnin:sim-1], probs = 0.975, na.rm = T)
estcoef[23] = quantile(holdomega2[burnin:sim-1], probs = 0.025, na.rm = T)
estcoef[24] = quantile(holdomega2[burnin:sim-1], probs = 0.975, na.rm = T)
estcoef[25] = as.numeric(estcoef[9]<holdmu[2,1]&estcoef[10]>holdmu[2,1])
estcoef[26] = as.numeric(estcoef[11]<holdmu[3,1]&estcoef[12]>holdmu[3,1])
estcoef[27] = as.numeric(estcoef[13]< holdomega1[1] &estcoef[14]> holdomega1[1])
estcoef[28] = as.numeric(estcoef[15]<holdmu[1,1]&estcoef[16]>holdmu[1,1])
estcoef[29] = as.numeric(estcoef[17]<holdmu[1,1]&estcoef[18]>holdmu[1,1])
#estcoef[30] = as.numeric(estcoef[19]<1&estcoef[20]>1)
estcoef[31] = as.numeric(estcoef[21]< holdpsi2[1] &estcoef[22]> holdpsi2[1])
estcoef[32] = as.numeric(estcoef[23]< holdomega2[1]&estcoef[24]>holdomega2[1])

## save results
fname = paste('params','coefsimVaccine.n',n,array_id2,'.txt',sep="")
write.table(params, file=fname, sep="\t", row.names=F, col.names=T)
fname2 = paste('PS','coefsimVaccine.n',n,array_id2,'.txt',sep="")
write.table(PS, file=fname2, sep="\t", row.names=F, col.names=T)
fname3 = paste('prentice','coefsim4Vaccine.n',n,array_id2,'.txt',sep="")
###write.table(prentice, file=fname3, sep="\t", row.names=F, col.names=T)
fname4 = paste('postpren','coefsimVaccine.n',n, array_id2,'.txt',sep="")
###write.table(pren.post, file=fname4, sep="\t", row.names=F, col.names=T)
fname5 = paste('naivemodels','coefsimVaccine.n',n,array_id2,'.txt',sep="")
###write.table(naiveresults, file=fname5, sep="\t", row.names=F, col.names=T)
fname6 = paste('estimatedcoef','coefsimVaccine.n',n,array_id2,'.txt',sep="")
write.table(estcoef, file=fname6, sep="\t", row.names=F, col.names=T)
fname7 = paste('covs','coefSimVaccine.n',n, array_id2,'.txt',sep="")
write.table(covs, file=fname7, sep="\t", col.names=T)
