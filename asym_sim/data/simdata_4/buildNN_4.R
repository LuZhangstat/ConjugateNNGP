setwd("H:/research/proj1")
setwd("./asym_sim/data/simdata_4")
rm(list = ls())
source("../NNMatrix.R")
load("simdata_4.RData")

## Build nearest neighbor matrix from spNNGP ##
M = 5 #10
library(spNNGP)
sigma.sq.IG <- c(as, bs)
cov.model <- "exponential"
theta.alpha <- c((ap + bp) / 2, bt / bs)
names(theta.alpha) <- c("phi", "alpha")

m.c <- spConjNNGP(Y ~ X[, 2], 
                  coords= coords, 
                  n.neighbors = M,
                  theta.alpha = theta.alpha, k.fold = 1,
                  n.omp.threads = 2, return.neighbors = T,
                  sigma.sq.IG = sigma.sq.IG, cov.model = cov.model)

NN.matrix <- NNMatrix(N, m.c$coords.ord, m.c$n.indx[-1])

par(mfrow = c(1, 1))
for(i in 200:220){
  plot(m.c$coords.ord) # all points
  points(m.c$coords.ord[1:i, , drop = FALSE], col = "grey", pch = 19) 
  
  # neighbors
  if (i < M) {dim = i} else {dim = M}
  for (j in 1:dim){
    points(m.c$coords.ord[NN.matrix$nearind[i - 1, j], , drop = FALSE], 
           col = "orange",  pch = 19)
  }
  points(m.c$coords.ord[i, , drop = FALSE], col = "blue", pch = 19) # the target point
  
  readline(prompt = "Pause. Press <Enter> to continue...")
}

save(list = ls(all.names = TRUE), file = "nngp_5_4.RData", envir = .GlobalEnv)

#save(list = ls(all.names = TRUE), file = "nngp_10_4.RData", envir = .GlobalEnv)







