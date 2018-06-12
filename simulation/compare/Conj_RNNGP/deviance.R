setwd("/Users/luzhang/Documents/github/ConjugateNNGP")
setwd("./simulation")
rm(list = ls())
library("sp")
library("rstan")
library("Matrix")
source("./src/functions.R")
load("./data/buildNN/nngp_10.RData")
load("./results/Conj_RNNGP.RData")

phi_cv <- 17.64706; deltasq_cv <- 0.09000

# recover a and b for sigma.sq
a <- m.c$sigma.sq.hat^2 / m.c$sigma.sq.var + 2
b <- m.c$sigma.sq.hat * (a - 1)
# expectation of ln(sigma.sq):
E_l_sigma.sq <- log(b) - digamma(a)
E_inv_sigma.sq <- a / b

D_ord <- dist(NN.matrix$coords.ord)
Cov_true <- sigma.sq * exp(- phi * as.matrix(D_ord)) + 
  tau.sq * diag(N)
chol_Cov_true <- chol(Cov_true)
logdet_true <- 2 * sum(log(diag(chol_Cov_true)))
X.ord <- X[NN.matrix$ord, ]

ind_x <-c(c(rep(2:M, times = 1:(M - 1)), rep(((M + 1) : N), each = M)), 1:N)
ind_y <- c(c(t(NN.matrix$NN_ind))[which(c(t(NN.matrix$NN_ind)) > 0)], 1:N)
phi_cv <- 17.64706; deltasq_cv <- 0.09000
AD <- getADstan2(neardist = NN.matrix$NN_dist,  neardistM = NN.matrix$NN_distM, 
                N = N, M = M, phi = phi_cv, deltasq = deltasq_cv) 
D <- AD[M + 1, ]
I_L <- sparseMatrix(i = ind_x, j = ind_y, dims = c(N, N), 
                    x = c(-na.omit(as.vector(AD[-(M + 1), ])), rep(1, N)))
inv_Cov_sim <- (t(I_L / D) %*% I_L) 
inv_Cps_Cq <-  inv_Cov_sim %*% Cov_true

inv_Cov_sim_XVBXT <- inv_Cov_sim %*% X.ord %*% m.c$beta.var %*% t(X.ord)

xb_xbs <- X.ord %*% (B - t(m.c$beta.hat))
chol_xb_xbs <- (I_L %*% xb_xbs) / sqrt(D)

KL_D <- sum(diag(inv_Cps_Cq)) * E_inv_sigma.sq + N * E_l_sigma.sq + 
  sum(log(D)) - logdet_true +  E_inv_sigma.sq * 
  (sum(diag(inv_Cov_sim_XVBXT)) + sum(chol_xb_xbs^2)) - N
KL_D <- 0.5 * KL_D
round(KL_D, 2)
#[1] 3.54
  