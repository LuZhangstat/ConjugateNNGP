setwd("/Users/luzhang/Documents/github/ConjugateNNGP") # set to the path of ConjugateNNGP
setwd("./SST_study")
rm(list = ls())
source("src/NNMatrix.R")
load("data/data/realdata.RData")


## Build nearest neighbor matrix from spNNGP ##
M = 10
library(spNNGP)
sigma.sq.IG <- c(as, bs)
theta.alpha <- c((ap + bp) / 2, bt / bs)
names(theta.alpha) <- c("phi", "alpha")

NN.matrix <- NNMatrix(coords = coords, n.neighbors = M,
                      n.omp.threads = 2, search.type = "brute")

save(list = ls(all.names = TRUE), file = "data/buildNN/nngp_10.RData", envir = .GlobalEnv)













