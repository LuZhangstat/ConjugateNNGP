# RMSPE 
setwd("/Users/luzhang/Documents/github/ConjugateNNGP")
setwd("./simulation")
rm(list = ls())
library("sp")
library("rstan")
library("Matrix")
library("fields")
load("./data/sorted_x_order/nngp_10.RData")
load("./results/conj_model2.RData")


M = 10
set.seed(1)
# get posterior samples
n_run <- 500; L <- length(pos.p.sigmasq)
sampleind <- sample.int(L, n_run)

Y_pred <- Matrix(0, nrow = (N_total - N), ncol = n_run)

for (l in (1):(N_total - N)){
  cat(l, "th new sample: \n")
  NNind_l <- sort(order(
    rdist(coords_test[c(1,l),], coords[c(1:N), ])[-1,])[1:M])
  NNdist <- dist(rbind(coords_test[l,], coords[NNind_l, ]))
  NNdistv <- c(NNdist)[1:M]
  NNdistM <- as.matrix(NNdist)[-1, -1]
  
  i = 0
  for (s in sampleind){
    i = i + 1
    beta_sim <- pos.p.beta[,s]
    sigmasq_sim <- pos.p.sigmasq[s]
    tausq_sim <- 0.05*pos.p.sigmasq[s]
    phi_sim <- 30
    y_xb_sim <- Y[NNind_l] - X[NNind_l, ]%*%beta_sim
    
    NNcorM <- exp(- phi_sim * NNdistM ) + (tausq_sim / sigmasq_sim) * diag(M) 
    chol_NNcorM <- t(chol(NNcorM))
    NNcorV <- exp(- phi_sim * NNdistv)
    v <- forwardsolve(chol_NNcorM, NNcorV)
    
    
    var_l <- sigmasq_sim * (1 - sum(v^2))
    
    mean_l <- X_test[l, ] %*% beta_sim + t(backsolve(t(chol_NNcorM), v)) %*% y_xb_sim
    Y_pred[l, i] <- mean_l + rnorm(1)* sqrt(var_l)
  }
}

MSPE <- mean((Y_test - rowMeans(Y_pred))^2)
RMSPE <- sqrt(MSPE); RMSPE

save(Y_pred, file = "compare/model2/Y_pred_MSPE.RData")


#----------------------- MSE of random effect --------------------#

set.seed(1)
n_run <- 500; L <- length(pos.p.sigmasq)
sampleind <- sample.int(L, n_run)

MSE.w <- mean(colSums((pos.p.w[, sampleind] - w[m.c$ord])^2)); MSE.w







