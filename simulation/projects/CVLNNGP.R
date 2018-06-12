setwd("/Users/luzhang/Documents/github/ConjugateNNGP") # set to the path of ConjugateNNGP
setwd("./simulation")
rm(list = ls())
library(spNNGP)
library(MBA)
library(fields)
library(MASS)
source("./src/NNMatrix.R")
source("./src/functions.R")
source("./src/util.R") # util.R in package spNNGP for rmspe and crps
library(foreach)
library(doParallel)
library(parallel)
load("./data/simdata/simdata.RData")
load("./data/buildNN/CVNNmatric.RData") 


## Parameters
sigma.sq.IG <- c(as, bs)
cov.model <- "exponential"
phi.deltasq <- as.matrix(expand.grid(3 / seq(0.05, 0.5, by = 0.02), 
                                       seq(0.01, 0.2, by = 0.02)))
colnames(phi.deltasq) <- c("phi", "deltasq")

g <- nrow(phi.deltasq)
k.fold.scores <- matrix(0, g, 1)
score.rule.names <- c("rmspe")
colnames(k.fold.scores) <- score.rule.names


## Data precess: (faster for building nearest neighbors)
#colnames(X) <- c("intercept", "projX", "projY")
ord <- order(coords[, 1]) # sort by longitude
coords <- coords[ord, ]
X <- X[ord, , drop = FALSE]
Y <- Y[ord]

set.seed(123)
## Cross-validation

t <- proc.time()
for (i in 1:k.fold){
  ho.id <- flds[[i]]
  X.ho <- X[ho.id, ]
  Y.ho <- Y[ho.id]
  coords.ho <- coords[ho.id, ]
  n.ho <- length(ho.id)
  
  X.mod <- X[-ho.id, ]
  Y.mod <- Y[-ho.id]
  coords.mod <- coords[-ho.id, ]
  n.mod <- nrow(X.mod)
  
  # index of nonzero element
  ind_x <-c(c(rep(2:M, times = 1:(M - 1)),
              rep(((M + 1):n.mod), each = M)), 1:n.mod)
  ind_y <- c(c(t(NN.matrix[[i]]$NN_ind))[
    which(c(t(NN.matrix[[i]]$NN_ind)) > 0)], 1:n.mod)
  
  ind_x_X <- rep(1:n.mod, P); ind_y_X <- rep(1:P, each = n.mod)
  ind_x_X_up <- c(ind_x_X, 1:n.mod)
  ind_y_X_up <- c(ind_y_X, (P + 1):(n.mod + P))
  
  ind_Au_X <- rep(1:n.ho, each = M)
  ind_Au_Y <- c(t(nn.mod.ho[[i]]$nn.idx))
  
  no_cores <- detectCores() - 2
  cl <- makeCluster(no_cores, type = "FORK", outfile = "./results/cvNNGP_debug.txt")
  clusterSetRNGStream(cl = cl, iseed = 1234)
  
  CVrmspe <- parSapply(cl, 1:g, function(j){
#  for(j in 1:g){
    phi <- phi.deltasq[j, 1]; deltasq <- phi.deltasq[j, 2]
    delta <- sqrt(deltasq)
    
    ## obtain A and D using C and N(i)
    AD <- getADstan(neardist = NN.matrix[[i]]$NN_dist,  
                    neardistM = NN.matrix[[i]]$NN_distM, 
                    N = n.mod, M = M, phi = phi) 
    D <- AD[M + 1, ]
    
    ## generate sparse matrix X* and Y* on S[-k]
    X_star_up <- 
      sparseMatrix(ind_x_X_up, ind_y_X_up,
                   x = c(as.vector(X.mod) / delta,
                         rep(1 / delta, n.mod)))
    X_star_down <- 
      sparseMatrix(i = ind_x, j = (ind_y + P), 
                   x = c(-na.omit(as.vector(AD[-(M + 1), ])), rep(1, n.mod))) / 
      sqrt(D)  
    
    #Dsqrtinv <- sparseMatrix(i = 1:N, j = 1:N, x = 1 / sqrt(D))
    X_star <- rbind(X_star_up, X_star_down)
    Y_star <- c(Y.mod / delta, rep(0, n.mod))
    
    ## get gamma_hat_k = (X*^T X*)^-1 X*^T y* by conjugate gradient ##
    X_Tstar_Y_star <- crossprod(X_star, Y_star)
    XTX_star <- crossprod(X_star, X_star)
    gamma_hat <- cgsparse(XTX_star, X_Tstar_Y_star[1:(n.mod + P)])
    
    ### sample prediction over holdout locations
    ## Build matrix A_u
    AD_ho <- getADstan_ho(neardist = nn.mod.ho[[i]]$nn.dists,  
                          neardistM = nn.mod.ho[[i]]$NN_distM, 
                          N = n.ho, M = M, phi = phi, deltasq = deltasq) 
    XU_Au <- cbind(X.ho, sparseMatrix(ind_Au_X, ind_Au_Y,
                                      x = as.vector(AD_ho[-(M + 1), ])))
    if(dim(XU_Au)[2] < length(gamma_hat)){
        gamma_hat <- gamma_hat[1:dim(XU_Au)[2]]
    }
    Y.ho.m <- XU_Au %*% gamma_hat
    return(c(rmspe(Y.ho, Y.ho.m@x), as.integer(j)))
  })
  k.fold.scores[CVrmspe[2, ], "rmspe"] <- 
    k.fold.scores[CVrmspe[2, ], "rmspe"] + CVrmspe[1, ]
#    k.fold.scores[j, "rmspe"] <- k.fold.scores[j, "rmspe"] + 
#      rmspe(Y.ho, Y.ho.m@x)
#  }
  stopCluster(cl)
  cat("\n", "--------------------", i/k.fold*100, "% done --------------------")
}
proc.time() - t
#user  system elapsed 
#0.237   0.095  12.186 
pick.ind <- which.min(k.fold.scores[, "rmspe"])
phi.deltasq[pick.ind, ]
#phi  deltasq 
#17.64706  0.09000 
save(file = "./results/CVscores.RData",
     list = c("k.fold.scores", "phi.deltasq"))

