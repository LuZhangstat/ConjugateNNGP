
setwd("./simulation")
rm(list = ls())
library("sp")
library("rstan")
library("Matrix")
load("./data/sorted_x_order/nngp_10.RData")
load("./results/response_GP")

Cov_true <- sigma.sq * exp(- phi * as.matrix(D_lit)) + tau.sq * diag(N)
chol_Cov_true <- chol(Cov_true)

# get posterior samples
post_sim <- extract(samples_full)

set.seed(1)
n_run <- 500; L <- length(post_sim$sigmasq)
sampleind <- sample.int(L, n_run)

t <- proc.time()
# calculate score and deviance
score <- rep(0, n_run); deviance <- rep(0, n_run)
simind = 0
for (s in sampleind){
  simind = simind + 1
  cat(s, "th iteration: \n")
  beta_sim <- post_sim$beta[s,]
  sigmasq_sim <- post_sim$sigmasq[s]
  tausq_sim <- post_sim$tausq[s]
  phi_sim <- post_sim$phi[s]
  
  Cov_sim <- sigmasq_sim * exp(- phi_sim * as.matrix(D_lit)) + 
    tausq_sim * diag(N)
  chol_Cov_sim <- chol(Cov_sim)
  y_xb_sim <- Y - X%*%beta_sim
  xb_xbs <- X%*%(B - beta_sim)
  chol_yxb_sim <- forwardsolve(t(chol_Cov_sim), y_xb_sim)
  chol_xb_xbs <- forwardsolve(t(chol_Cov_sim), xb_xbs)
  
  inv_Cps_Cq <-  chol2inv(chol_Cov_sim)%*% Cov_true
  
  score[simind] <- - 2 * sum(log(diag(chol_Cov_sim))) - sum(chol_yxb_sim^2)
  
  deviance[simind] <- sum(diag(inv_Cps_Cq)) -
    determinant(inv_Cps_Cq, logarithm = TRUE)$modulus +
    sum(chol_xb_xbs^2)- N
  
}
proc.time() - t

save(score, deviance, 
     file = "compare/response_GP/score_dev.RData")
