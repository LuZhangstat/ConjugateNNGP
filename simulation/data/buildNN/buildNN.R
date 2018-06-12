setwd(".")
rm(list = ls())
source("./simulation/src/NNMatrix.R")
load("./simulation/data/simdata/simdata.RData")

## Build nearest neighbor matrix from spNNGP ##
M = 10

NN.matrix <- NNMatrix(coords = coords, n.neighbors = M,
                      n.omp.threads = 2, search.type = "tree")

par(mfrow = c(1, 1))
Check_Neighbors(NN.matrix$coords.ord, n.neighbors = M, NN.matrix, 600)

save(list = ls(all.names = TRUE), 
     file = "./simulation/data/buildNN/nngp_10.RData", envir = .GlobalEnv)







