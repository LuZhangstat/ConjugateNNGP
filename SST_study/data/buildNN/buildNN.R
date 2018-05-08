setwd("/Users/luzhang/Documents/github/ConjugateNNGP") # set to the path of ConjugateNNGP
setwd("./SST_study")
rm(list = ls())
source("projects/NNMatrix.R")
load("data/data/realdata.RData")


## Build nearest neighbor matrix from spNNGP ##
M = 10
library(spNNGP)
sigma.sq.IG <- c(as, bs)
cov.model <- "exponential"
theta.alpha <- c((ap + bp) / 2, bt / bs)
names(theta.alpha) <- c("phi", "alpha")

NN.matrix <- NNMatrix(coords = coords, n.neighbors = M,
                      n.omp.threads = 2, search.type = "brute")

save(list = ls(all.names = TRUE), file = "data/buildNN/nngp_10.RData", envir = .GlobalEnv)


#par(mfrow = c(1, 1))
#for(i in 200:220){
#  plot(m.c$coords.ord[1:2000,]) # all points
#  points(m.c$coords.ord[1:i, , drop = FALSE], col = "grey", pch = 19)
  
  # neighbors
#  if (i < M) {dim = i} else {dim = M}
#  for (j in 1:dim){
#    points(m.c$coords.ord[NN.matrix$nearind[i - 1, j], , drop = FALSE],
#           col = "orange",  pch = 19)
#  }
#  points(m.c$coords.ord[i, , drop = FALSE], col = "blue", pch = 19) # the target point
  
#  readline(prompt = "Pause. Press <Enter> to continue...")
#}










