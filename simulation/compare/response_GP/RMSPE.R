# RMSPE
setwd("./simulation")
rm(list = ls())
library("sp")
library("rstan")
library("Matrix")
library("fields")
load("./results/response_GP")
load("./data/simdata/simdata.RData")

M = 10
set.seed(1)

# get posterior samples
post_sim <- extract(samples_full)
n_run <- 500; L <- length(post_sim$sigmasq)
sampleind <- sample.int(L, n_run)

Y_pred <- Matrix(0, nrow = (N_total - N), ncol = n_run)
i = 0
for (s in sampleind){
  i = i + 1
  cat(s, "th posterior simulation: \n")
  beta_sim <- post_sim$beta[s,]
  sigmasq_sim <- post_sim$sigmasq[s]
  tausq_sim <- post_sim$tausq[s]
  phi_sim <- post_sim$phi[s]
  y_xb_sim <- Y - X[1:N, ]%*%beta_sim
    
  corM <- exp(- phi_sim * D[1:N, 1:N]) + (tausq_sim / sigmasq_sim) * diag(N) 
  chol_corM <- t(chol(corM))
    
  for (l in (1):(N_total - N)){
    distv <- rdist(coords_test[c(1, l), ], coords)[-1, ]
    corV <- exp(- phi_sim * distv)
    v <- forwardsolve(chol_corM, corV)
    
    var_l <- (sigmasq_sim + tausq_sim) - sigmasq_sim * sum(v^2)
    
    mean_l <- X_test[l, ]%*%beta_sim + t(backsolve(t(chol_corM), v)) %*% y_xb_sim
    Y_pred[l, i] <- mean_l + rnorm(1)* sqrt(var_l)
  }
}

MSPE <- mean((Y_test - rowMeans(Y_pred))^2)
RMSPE <- sqrt(MSPE); RMSPE

save(Y_pred, file = "compare/response_GP/Y_pred_MSPE.RData")




