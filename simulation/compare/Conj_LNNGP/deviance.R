setwd("/Users/luzhang/Documents/github/ConjugateNNGP")
setwd("./simulation")
rm(list = ls())
library("sp")
library("rstan")
library("Matrix")
source("./src/functions.R")
load("./data/buildNN/nngp_10.RData")
load("./results/Conj_LNNGP_paral.RData")

# get posterior samples
set.seed(1)
n_run <- 300; L <- length(Conj_LNNGP_paral)
sampleind <- sample.int(L, n_run)


D_ord <- dist(NN.matrix$coords.ord)
Cov_true <- sigma.sq * exp(- phi * as.matrix(D_ord)) + 
  tau.sq * diag(N)
chol_Cov_true <- chol(Cov_true)
logdet_true <- 2 * sum(log(diag(chol_Cov_true)))
X.ord <- X[NN.matrix$ord, ]

ind_x <-c(c(rep(2:M, times = 1:(M - 1)), rep(((M + 1) : N), each = M)), 1:N)
ind_y <- c(c(t(NN.matrix$NN_ind))[which(c(t(NN.matrix$NN_ind)) > 0)], 1:N)
phi_cv <- 17.64706; deltasq_cv <- 0.09000
AD <- getADstan(neardist = NN.matrix$NN_dist,  neardistM = NN.matrix$NN_distM, 
                N = N, M = M, phi = phi_cv) 
D <- AD[M + 1, ]
I_L <- sparseMatrix(i = ind_x, j = ind_y, dims = c(N, N), 
                    x = c(-na.omit(as.vector(AD[-(M + 1), ])), rep(1, N)))
inv_Cov_sim <- (t(I_L / D) %*% I_L)

t <- proc.time()
# calculate score and deviance
KL_D <- rep(0, n_run)
simind = 0
for (s in sampleind){
  simind = simind + 1
   
  sigmasq_sim <- Conj_LNNGP_paral[[s]]$pos.p.sigmasq
  tausq_sim <- sigmasq_sim * deltasq_cv
  beta_sim <- Conj_LNNGP_paral[[s]]$pos.p.beta
  
  Cov_sim <- sigmasq_sim *chol2inv(chol(inv_Cov_sim))
  diag(Cov_sim) = diag(Cov_sim) + tausq_sim
  chol_Cov_sim <- chol(Cov_sim)
  xb_xbs <- X.ord %*% (B - beta_sim)
  chol_xb_xbs <- forwardsolve(t(chol_Cov_sim), xb_xbs)
  inv_Cps_Cq <-  chol2inv(chol_Cov_sim)%*% Cov_true
  
  KL_D[simind] <- 0.5 * sum(diag(inv_Cps_Cq)) + 
    sum(log(diag(chol_Cov_sim))) - 0.5 * logdet_true + 
    sum(chol_xb_xbs^2) / 2 - N / 2
  cat(simind, "th iteration: cost", (proc.time() - t)[3], "\n")
}
proc.time() - t

round(mean(KL_D), 2)
#[1] 3.58
round(quantile(KL_D, c(0.5, 0.025, 0.975)), 2)
#50%  2.5% 97.5% 
#3.03  1.27  8.56 

#----------------------- MSE(w) = mean[(w_fit - w)^2]--------------------------#
simind = 0; MSE_w = rep(0, n_run)
w.ord <- w[NN.matrix$ord]
for (s in sampleind){
  simind = simind + 1
  MSE_w[simind] <- sum((Conj_LNNGP_paral[[s]]$pos.p.w - w.ord)^2)
}
round(mean(MSE_w), 2)
#[1] 313.28
round(quantile(MSE_w, c(0.5, 0.025, 0.975)), 2)
#50%   2.5%  97.5% 
#292.20 258.96 483.75
