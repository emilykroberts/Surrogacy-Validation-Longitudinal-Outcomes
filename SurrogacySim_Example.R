# load required packages
library(corpcor); library(bayesSurv); library(MASS); library(GenKern); library(coda); library(psych)
library(mvtnorm); library(rootSolve); library(Matrix); library(MCMCpack); library(LaplacesDemon)
library(reshape2); library(tidyr); library(ggplot2); library(pracma); library(lme4); library(lcmm)
library(HardyWeinberg); library(ExtDist)

setwd("./R")
# keep source files in a folder called R directory
for(i in 1:(length(list.files(pattern = "\\.R$")))){
    source(list.files(pattern = "\\.R$")[i]
 )
}



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
misspecify = F # lcmix pacakage and mv gamma no longer available
r=3

# generating values
alpha1 = 2; psi2 = 0
beta0 = 2; omega1 = 0
beta1 = 3.1; omega2 = 0
eps1 = 1; eps2 = 1; eps3 = 1; tau4 = 0.5; eta4 = 1
theta10 = 0.15; theta11 = 0.7; thetaT = 0.2142857
sigmae = 0.3; tau = sigmae^2 # represents ei
R = matrix(rep(1,3*3),3,3); R[1,2] = R[2,1] = theta10; R[1,3] = R[3,1] = theta11; R[2,3] = R[3,2] = thetaT

mu = c(alpha, beta0, beta1)
### generate longitudinal T0 and T1 = mu + ui + eij #ui shared across outcomes per person

samp = generate_data(n = n, theta10 = theta10, theta11 = theta11, thetaT = thetaT, alpha1 = alpha1, 
psi2 = psi2,
                         beta0 = beta0, omega1 = omega1,
                         beta1 = beta1, omega2 = omega2,
                         eps1 = eps1, eps2 = eps2, eps3 = eps3, 
                         sigmae = sigmae, misspecify = misspecify, r = 3)
                         
# needed hyperparameters
a = b = 0.1

# set up matrices to save results
# save data corresponding to this simulation on cluster

trt = c(rep(0,n)); trt[(1*n-(n/2-1)):(n*1)] = 1
ST = samp[[1]]

uitruesamp = samp[[2]]

holdmu = matrix(rep(0,3*SIM),3,SIM); holdS = array(rep(0,3*3*SIM),dim=c(3,3,SIM))
holdR = array(rep(0,3*3*SIM),dim = c(3,3,SIM)); holdR[,,1] = diag(c(1,1,1)); holdS[,,1] = diag(c(1,1,1))
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
bivec = rep(0, n)

x_1 = z_1 = matrix(c(c(1, 1, 1, 1)), ncol = 2, byrow = F)
Zi = Zi0 = Zi1 = rep(list(z_1), n)
X1 = rep(list(x_1, n))

for(i in 1:n){
  
  X1[[i]] = matrix(c(c(1, 1, 1, 1, 1, 1)), ncol = 2, byrow = T)
  Zi1[[i]] = matrix(c(c(1, 1, 1)), ncol = 3, byrow = T)
  Zi[[i]] = matrix(c(c(1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1)), ncol = 2, byrow = T)
  
}

ST = cbind(ST, rep(0, n), rep(0, n))

# initial values
r13init = 0.5
r12init = 0.1
r23init = 0.1

# start simulation

run_mcmc(samp = ST, Xi = rep(0,n), n = n, SIM = SIM, r = r, condindfit = condindfit,
                    prior = prior, r13init = r13init, r12init = r12init, r23init = r23init,
                    b0init = rep(0, n/2), b1init = rep(0, n/2))
  
  