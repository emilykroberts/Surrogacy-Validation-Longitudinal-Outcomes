
#' Write results
#'
#' @description Write results from mcmc
#'
#' @param params parameters from mcmc
#' @param write logical value of saving results
#'
#' @return results
#'
#' @examples
#' example(final_results(params, write = F))
#'
final_results = function(params, write){
  ## save results
  param = params$params

  n = params$args$n
  burnin = params$args$burnin
  sim = params$args$SIM


  PS = data.frame(beta0SE = numeric(1), beta1SE = numeric(1), alpha1SE = numeric(1),
                  psi1SE = numeric(1), omega1SE = numeric(1), psi2SE = numeric(1), omega2SE = numeric(1),
                  sigs0 = numeric(1), SEsigs0 = numeric(1), sigs1 = numeric(1), SEsigs1 = numeric(1),
                  sigmat0 = numeric(1), SEsigmat0 = numeric(1),
                  ps = numeric(1), SEps = numeric(1),
                  p10 = numeric(1),SEp10 = numeric(1),p11 = numeric(1),
                  SEp11 = numeric(1),S0l = numeric(1), S0u = numeric(1),
                  S1l = numeric(1), S1u = numeric(1), T0l = numeric(1),T0u = numeric(1),
                  psl = numeric(1), PSu = numeric(1),
                  p10l = numeric(1),p10u = numeric(1),
                  p11l = numeric(1),p11u = numeric(1),

                  b0 = numeric(1), b0SE = numeric(1), b1 = numeric(1), b1SE = numeric(1),
                  beta0mean = numeric(1), beta1mean = numeric(1), omega1mean = numeric(1),
                  alpha1mean = numeric(1),psi1mean = numeric(1),psi2mean = numeric(1),omega2mean = numeric(1),
                  beta0l = numeric(1),beta0u = numeric(1), beta1l = numeric(1),beta1u = numeric(1),omega1l = numeric(1),
                  omega1u = numeric(1),
                  alpha1l = numeric(1),alpha1u = numeric(1),psi1l = numeric(1),psi1u = numeric(1),
                  psi2l = numeric(1),psi2u = numeric(1),omega2l = numeric(1),omega2u = numeric(1),
                  mean_int = numeric(1), SEmean_int = numeric(1),L_int = numeric(1), U_int = numeric(1),mean_sl = numeric(1),
                  SEmean_sl = numeric(1),L_sl = numeric(1),U_sl = numeric(1), gamma0g0 = numeric(1), gamma1g0 = numeric(1))



  PS[1] = sqrt(var(param$beta0[burnin:(sim-1)]))
  PS[2] = sqrt(var(param$beta1[burnin:(sim-1)]))
  PS[3] = sqrt(var(param$alpha1[burnin:(sim-1)]))
  PS[4] = sqrt(var(param$psi1[burnin:(sim-1)]))
  PS[5] = sqrt(var(param$omega1[burnin:(sim-1)]))
  PS[6] = sqrt(var(param$psi2[burnin:(sim-1)]))
  PS[7] = sqrt(var(param$omega2[burnin:(sim-1)]))

  PS[8] = mean(param$s1[burnin:(sim-1)])
  PS[9] = sqrt(var(param$s1[burnin:(sim-1)]))
  PS[10] = mean(param$s2[burnin:(sim-1)])
  PS[11] = sqrt(var(param$s2[burnin:(sim-1)]))
  PS[12] = mean(param$s3[burnin:(sim-1)])
  PS[13] = sqrt(var(param$s3[burnin:(sim-1)]))

  PS[14] = mean(param$r12[burnin:(sim-1)])
  PS[15] = sqrt(var(param$r12[burnin:(sim-1)]))
  PS[16] = mean(param$r13[burnin:(sim-1)])
  PS[17] = sqrt(var(param$r13[burnin:(sim-1)]))
  PS[18] = mean(param$r23[burnin:(sim-1)])
  PS[19] = sqrt(var(param$r23[burnin:(sim-1)]))

  PS[20] = quantile(param$s1[burnin:(sim-1)], probs = 0.025, na.rm = TRUE)
  PS[21] = quantile(param$s1[burnin:(sim-1)], probs = 0.975, na.rm = TRUE)
  PS[22] = quantile(param$s2[burnin:(sim-1)], probs = 0.025, na.rm = TRUE)
  PS[23] = quantile(param$s2[burnin:(sim-1)], probs = 0.975, na.rm = TRUE)
  PS[24] = quantile(param$s3[burnin:(sim-1)], probs = 0.025, na.rm = TRUE)
  PS[25] = quantile(param$s3[burnin:(sim-1)], probs = 0.975, na.rm = TRUE)

  PS[26] = quantile(param$r12[burnin:(sim-1)], probs = 0.025, na.rm = TRUE)
  PS[27] = quantile(param$r12[burnin:(sim-1)], probs = 0.975, na.rm = TRUE)
  PS[28] = quantile(param$r13[burnin:(sim-1)], probs = 0.025, na.rm = TRUE)
  PS[29] = quantile(param$r13[burnin:(sim-1)], probs = 0.975, na.rm = TRUE)
  PS[30] = quantile(param$r23[burnin:(sim-1)], probs = 0.025, na.rm = TRUE)
  PS[31] = quantile(param$r23[burnin:(sim-1)], probs = 0.975, na.rm = TRUE)

  PS[32] = mean(param$b0[burnin:(sim-1)])
  PS[33] = sqrt(var(param$b0[burnin:(sim-1)]))
  PS[34] = mean(param$b1[burnin:(sim-1)])
  PS[35] = sqrt(var(param$b1[burnin:(sim-1)]))
  PS[36] = mean(param$beta0[burnin:sim-1], na.rm = T)
  PS[37] = mean(param$beta1[burnin:sim-1], na.rm = T)
  PS[38] = mean(param$omega1[burnin:sim-1], na.rm = T)
  PS[39] = mean(param$alpha1[burnin:sim-1], na.rm = T)
  PS[40] = mean(param$psi1[burnin:sim-1], na.rm = T)
  PS[41] = mean(param$psi2[burnin:sim-1], na.rm = T)
  PS[42] = mean(param$omega2[burnin:sim-1], na.rm = T)
  PS[43] = quantile(param$beta0[burnin:sim-1], probs = 0.025, na.rm = T)
  PS[44] = quantile(param$beta0[burnin:sim-1], probs = 0.975, na.rm = T)
  PS[45] = quantile(param$beta1[burnin:sim-1], probs = 0.025, na.rm = T)
  PS[46] = quantile(param$beta1[burnin:sim-1], probs = 0.975, na.rm = T)
  PS[47] = quantile(param$omega1[burnin:sim-1], probs = 0.025, na.rm = T)
  PS[48] = quantile(param$omega1[burnin:sim-1], probs = 0.975, na.rm = T)
  PS[49] = quantile(param$alpha1[burnin:sim-1], probs = 0.025, na.rm = T)
  PS[50] = quantile(param$alpha1[burnin:sim-1], probs = 0.975, na.rm = T)
  PS[51] = quantile(param$psi1[burnin:sim-1], probs = 0.025, na.rm = T)
  PS[52] = quantile(param$psi1[burnin:sim-1], probs = 0.975, na.rm = T)
  PS[53] = quantile(param$psi2[burnin:sim-1], probs = 0.025, na.rm = T)
  PS[54] = quantile(param$psi2[burnin:sim-1], probs = 0.975, na.rm = T)
  PS[55] = quantile(param$omega2[burnin:sim-1], probs = 0.025, na.rm = T)
  PS[56] = quantile(param$omega2[burnin:sim-1], probs = 0.975, na.rm = T)

  PS[57] = mean(param$int[burnin:sim-1], na.rm = T)
  PS[58] = sqrt(var(param$int[burnin:sim-1], na.rm = T))
  PS[59] = quantile(param$int[burnin:sim-1], probs = 0.025, na.rm = T)
  PS[60] = quantile(param$int[burnin:sim-1], probs = 0.975, na.rm = T)
  PS[61] = mean(param$slope[burnin:sim-1], na.rm = T)
  PS[62] = sqrt(var(param$slope[burnin:sim-1], na.rm = T))
  PS[63] = quantile(param$slope[burnin:sim-1], probs = 0.025, na.rm = T)
  PS[64] = quantile(param$slope[burnin:sim-1], probs = 0.975, na.rm = T)
  PS[65] = as.numeric(PS[59]<0 & PS[60]>=0)
  PS[66] = as.numeric(PS[63]>0)

  print(PS)
  ## save results
  if(write){
    fname2 = paste('PS','coefsimVaccine','.txt',sep="")
    write.table(PS, file=fname2, sep="\t", row.names=F, col.names=T)
  }

}
