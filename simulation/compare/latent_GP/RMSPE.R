# RMSPE 

setwd("./simulation")
rm(list = ls())
library("sp")
library("rstan")
library("Matrix")
library("fields")

load("./results/latent_GP")
load("./data/simdata/simdata.RData")

M = 10
set.seed(1)

# get posterior samples
post_sim <- extract(samples_w_f)
n_run <- 500; L <- length(post_sim$sigmasq)
set.seed(1)
sample_ind <- sample(L, n_run)

Y_pred <- Matrix(0, nrow = (N_total - N), ncol = n_run)
for (ind in 1:n_run){
  s <- sample_ind[ind]
  cat(ind, "th iteration, ", s, "th posterior simulation: \n")
  beta_sim <- post_sim$beta[s,]
  sigmasq_sim <- post_sim$sigmasq[s]
  tausq_sim <- post_sim$tausq[s]
  phi_sim <- post_sim$phi[s]
  w_sim <- post_sim$w[s,]
  
  w_b1_sim <- w_sim - beta_sim[1]
  
  corM <- exp(- phi_sim * D[1:N, 1:N]) 
  chol_corM <- t(chol(corM))
  
  for (l in (1):(N_total - N)){
    distv <- rdist(coords_test[c(1, l), ], coords)[-1, ]
    corV <- exp(- phi_sim * distv)
    v <- forwardsolve(chol_corM, corV)
    
    var_l <- sigmasq_sim * (1 - sum(v^2))
    mean_l <- beta_sim[1] + t(backsolve(t(chol_corM), v)) %*% w_b1_sim
    
    w_new <- mean_l + rnorm(1)* sqrt(var_l)
      
    Y_pred[l, ind] <- X_test[l, 2]*beta_sim[2] + w_new + rnorm(1)* sqrt(tausq_sim)
  }
}

MSPE <- mean((Y_test - rowMeans(Y_pred))^2)
RMSPE <- sqrt(MSPE); RMSPE

save(Y_pred, file = "compare/latent_GP/Y_pred_MSPE.RData")

#----------------------- MSE of random effect --------------------#
set.seed(1)
n_run <- 500; L <- length(post_sim$sigmasq)
sampleind <- sample.int(L, n_run)

MSE.w <- mean(colSums((t(post_sim$w[sampleind,]- post_sim$beta[sampleind,1]) - 
w)^2)); MSE.w





