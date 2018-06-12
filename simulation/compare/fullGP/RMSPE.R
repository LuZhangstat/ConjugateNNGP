# RMSPE 

setwd("./simulation")
rm(list = ls())
library("spBayes")
load("./data/buildNN/nngp_10.RData")
load("./results/fullGP.RData")

burn.in = 10101
set.seed(123)
p.s <- spPredict(m.fullGP, pred.covars = X_test, pred.coords = coords_test, 
                 start = burn.in, thin = 33)

MSPE <- mean((Y_test - rowMeans(p.s$p.y.predictive.samples))^2)
RMSPE <- sqrt(MSPE); round(RMSPE, 2)

