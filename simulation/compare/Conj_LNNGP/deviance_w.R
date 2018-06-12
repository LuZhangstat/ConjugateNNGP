
setwd("./simulation")
rm(list = ls())
library("sp")
library("rstan")
library("Matrix")
load("./data/buildNN/nngp_10.RData")
load("./results/Conj_LNNGP_paral.RData")

phi_cv <- 17.64706; deltasq_cv <- 0.09000
y.ord <- Y[NN.matrix$ord]
X.ord <- X[NN.matrix$ord, ]
XBw.ord <- X.ord %*% B + w[NN.matrix$ord]

# get posterior samples
set.seed(1)
n_run <- 300; L <- length(Conj_LNNGP_paral)
sampleind <- sample.int(L, n_run)

# calculate score and deviance
KL_D <- rep(0, n_run)
simind = 0
for (s in sampleind){
  simind = simind + 1
  #cat(s, "th iteration: \n")
  sigmasq_sim <- Conj_LNNGP_paral[[s]]$pos.p.sigmasq
  tausq_sim <- deltasq_cv * sigmasq_sim
  
  xbw_xbw_sim <- XBw.ord - X.ord %*% Conj_LNNGP_paral[[s]]$pos.p.beta - 
    Conj_LNNGP_paral[[s]]$pos.p.w
  
  # KL-D = 0.5(-Nlog(tau.sq) - N + log(tausq_sim))
  KL_D[simind] <- 0.5 * (-N - N * log(tau.sq) + N * log(tausq_sim) + N * 
    tau.sq / tausq_sim + sum(xbw_xbw_sim^2) / tausq_sim)
}


round(mean(KL_D), 2)
round(quantile(KL_D, c(0.5, 0.025, 0.975)), 2)

