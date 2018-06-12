
setwd("./simulation")
rm(list = ls())
library("sp")
library("rstan")
library("Matrix")
load("./data/buildNN/nngp_10.RData")
load("./results/fullGP.RData")

XBw <- X %*% B + w

# get posterior samples
set.seed(1)
n_run <- 300; L <- 10000
sampleind <- sample.int(L, n_run)
# calculate score and deviance
KL_D <- rep(0, n_run)
simind = 0
for (s in sampleind){
  simind = simind + 1
  #cat(s, "th iteration: \n")
  tausq_sim <- m.fullGP$p.theta.samples[s + 10000, "tau.sq"]
  
  xbw_xbw_sim <- XBw - X %*% m.fullGP$p.beta.recover.samples[s, ] - 
    m.fullGP$p.w.recover.samples[, s]
  
  # KL-D = 0.5(-Nlog(tau.sq) - N + log(tausq_sim))
  KL_D[simind] <- 0.5 * (-N - N * log(tau.sq) + N * log(tausq_sim) + N * 
                           tau.sq / tausq_sim + sum(xbw_xbw_sim^2) / tausq_sim)
}

round(mean(KL_D), 2)
round(quantile(KL_D, c(0.5, 0.025, 0.975)), 2)
