setwd("/Users/luzhang/Documents/github/ConjugateNNGP")
setwd("./simulation")
rm(list = ls())
library("sp")
library("rstan")
library("Matrix")
load("./data/buildNN/nngp_10.RData")
load("./results/fullGP.RData")

# get posterior samples
set.seed(1)
n_run <- 10000; L <- 10000
sampleind <- sample.int(L, n_run)

D_train <- dist(coords)
Cov_true <- sigma.sq * exp(- phi * as.matrix(D_train)) + 
  tau.sq * diag(N)
chol_Cov_true <- chol(Cov_true)
logdet_true <- 2 * sum(log(diag(chol_Cov_true)))


t <- proc.time()
# calculate score and deviance
KL_D <- rep(0, n_run)
simind = 0
for (s in sampleind){
  simind = simind + 1
  #cat(s, "th iteration: \n")
  tausq_sim <- m.fullGP$p.theta.samples[s + 10000, "tau.sq"]
  sigmasq_sim <- m.fullGP$p.theta.samples[s + 10000, "sigma.sq"]
  phi_sim <- m.fullGP$p.theta.samples[s + 10000, "phi"]
  beta_sim <- m.fullGP$p.beta.recover.samples[s, ]
  
  Cov_sim <- sigmasq_sim * exp(-phi_sim * as.matrix(D_train))
  diag(Cov_sim) <- diag(Cov_sim) + tausq_sim
  chol_Cov_sim <- chol(Cov_sim)
  xb_xbs <- X %*% (B - beta_sim)
  chol_xb_xbs <- forwardsolve(t(chol_Cov_sim), xb_xbs)
  inv_Cps_Cq <-  chol2inv(chol_Cov_sim) %*% Cov_true
  
  KL_D[simind] <- 0.5 * sum(diag(inv_Cps_Cq)) + 
    sum(log(diag(chol_Cov_sim))) - 0.5 * logdet_true + 
    sum(chol_xb_xbs^2) / 2 - N / 2
}
proc.time() - t

round(mean(KL_D), 2)
#[1] 4.45
round(quantile(KL_D, c(0.5, 0.025, 0.975)), 2)
#50%  2.5% 97.5% 
#3.93  1.16  9.95 


#----------------------- MSE(w) = mean[(w_fit - w)^2]--------------------------#
simind = 0; MSE_w = rep(0, n_run)
for (s in sampleind){
  simind = simind + 1
  MSE_w[simind] <- sum((m.fullGP$p.w.recover.samples[, s] - w)^2)
}
round(mean(MSE_w), 2)
#[1] 297.45
round(quantile(MSE_w, c(0.5, 0.025, 0.975)), 2)
#50%   2.5%  97.5% 
#284.63 231.62 444.79 

# compare to all posterior samples
#[1] 295.4

#50%   2.5%  97.5% 
#281.82 228.70 450.96


