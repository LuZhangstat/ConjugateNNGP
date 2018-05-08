setwd("/Users/luzhang/Documents/github/ConjugateNNGP") # set to the path of ConjugateNNGP
setwd("./SST_study")
rm(list = ls())
library(spNNGP)
library(MBA)
library(fields)
load("./data/buildNN/nngp_10.RData")
load("./data/data/CVNNmatric.RData")
source("./projects/NNMatrix.R")
source("./projects/functions.R")
source("./projects/util.R") # util.R in package spNNGP for rmspe and crps


## Parameters
sigma.sq.IG <- c(as, bs)
cov.model <- "exponential"
phi.deltasq <- as.matrix(expand.grid(length(seq(3, 21, by = 2)), 
                                     seq(0.001, 0.02, by = 0.001)))
colnames(phi.deltasq) <- c("phi", "deltasq")

k.fold = 5; g <- nrow(phi.deltasq)
k.fold.scores <- matrix(0, g, 2)
score.rule.names <- c("rmspe", "crps")
colnames(k.fold.scores) <- score.rule.names


## Data precess: (faster for building nearest neighbors)
colnames(X) <- c("intercept", "projX", "projY")
ord <- order(coords[, 1]) # sort by longitude
coords <- coords[ord, ]
X <- X[ord, , drop = FALSE]
Y <- Y[ord]


## Build NN matrixs for Cross Validation
set.seed(123)
require(caret)
flds <- createFolds(1:N, k = k.fold, list = TRUE, returnTrain = FALSE)

#NN.matrix <- list()
#nn.mod.ho <- list()

#for (i in 1: k.fold) {
#  cat("\n", "Building NN matrics for ", i, "th folder")
#  t <- proc.time()
#  ho.id <- flds[[i]]
#  coords.ho <- coords[ho.id, ]
#  coords.mod <- coords[-ho.id, ]
  
  # Build NN index for training data in ith folder
#  NN.matrix[[i]] <- NNMatrix(coords = coords.mod, n.neighbors = M, 
#                             n.omp.threads = 2, search.type = "brute")
  
  # Build nearest neighbor index for holdout locations
#  nn.mod.ho[[i]] <- nn2(coords.mod, coords.ho, k = M)
#  nn.mod.ho[[i]]$NN_distM <- 
#    sapply(1:nrow(coords.ho), ho_dist, nn.mod.ho[[i]]$nn.idx, coords.mod)
#  cat("\t", " takes time: ", proc.time() - t)
#}

#save(file = "./data/data/CVNNmatric.RData", 
#     list = c("flds", "NN.matrix", "nn.mod.ho"))

set.seed(123)
## Cross-validation 
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
  
  for (j in 1:g){
    t <- proc.time()
    phi = phi.deltasq[j, 1]; deltasq <- phi.deltasq[j, 2]
    # index of nonzero element in 
    ind_x <-c(c(rep(2:M, times = 1:(M - 1)), 
                    rep(((M + 1) : n.mod), each = M)), 1:n.mod)
    ind_y <- c(c(t(NN.matrix[[i]]$NN_ind))[
      which(c(t(NN.matrix[[i]]$NN_ind)) > 0)], 1:n.mod)

    ## obtain A and D using C and N(i)
    AD <- getADstan(neardist = NN.matrix[[i]]$NN_dist,  
                    neardistM = NN.matrix[[i]]$NN_distM, 
                    N = n.mod, M = M, phi = phi) 
    D <- AD[, M + 1]
    
    ## generate sparse matrix X** and Y** 
    ind_x_X <- rep(1:n.mod, P); ind_y_X <- rep(1:P, each = n.mod)
    ind_x_X_up <- c(ind_x_X, 1:n.mod)
    ind_y_X_up <- c(ind_y_X, (P + 1):(n.mod + P))
    
    X_star_star_up <- 
      sparseMatrix(ind_x_X_up, ind_y_X_up,
                   x = c(1 / sqrt(deltasq) * c(X.mod),
                         rep(1 / sqrt(deltasq), n.mod)))
    X_star_star_down1 <- 
      sparseMatrix(i = ind_x, j = (ind_y + P), 
                   x = c( - c(t(AD[, 1:M]))[
                     which(!is.na(t(AD[, 1:M])))], rep(1, n.mod)))
    
    #Dsqrtinv <- sparseMatrix(i = 1:N, j = 1:N, x = 1 / sqrt(D))
    X_star_star_down <- X_star_star_down1 / sqrt(D)
    X_star_star <- rbind(X_star_star_up, X_star_star_down)
    Y_star_star <- c(Y.mod / sqrt(deltasq), rep(0, n.mod))
    
    ## get theta = (X**^T X**)^-1 X**^T y** by conjugate gradient ##
    
    X_Tstar_Y_star_star <- t(X_star_star) %*% Y_star_star
    XTX_star_star <- t(X_star_star) %*% X_star_star
    theta_hat <- cgsparse(XTX_star_star[1: (n.mod + P), 1:(n.mod + P)],
                          X_Tstar_Y_star_star[1:(n.mod + P)])
    
    ## obtain a* and b* by equation for the posterior IG(a*, b*) for sigmasq ##
    b_star_star <- bs + 0.5 * sum((Y_star_star - (X_star_star %*% theta_hat))^2) 
    a_star_star <- as + 0.5 * n.mod
    
    ## sample sigmasq from IG(a*, b*) ##
    sigmasq <- 1 / rgamma(1, shape = a_star_star, rate = b_star_star)
    cat("\t", sigmasq, "\t")
    
    #### sample theta given sigmasq, \hat{\theta} and X** ####
    ##i. sample u ~ N(0, sigmasq * I_{2n+p})$ ##
    
    u <- rnorm(2 * n.mod, 0) * sqrt(sigmasq)
    
    ##ii. get v = (X**T X**)^{-1} X**^Tu$ by conjugate gradient
    
    X_star_u <- t(X_star_star) %*% u
    v <- cgsparse(XTX_star_star[1: (n.mod + P), 1:(n.mod + P)],
                  X_star_u[1:(n.mod + P)])
    
    ##iii. get posterior sample theta = \hat{\theta} + v ##
    sample_theta <- theta_hat + v
    
    ### sample prediction over holdout locations
    ## Build matrix C and D
    ind_C.ho_X <- rep(1:n.ho, each = M)
    ind_C.ho_Y <- c(t(nn.mod.ho[[i]]$nn.idx))
    
    ## obtain A and D using C and N(i)
    AD_ho <- getADstan_ho(neardist = nn.mod.ho[[i]]$nn.dists,  
                    neardistM = nn.mod.ho[[i]]$NN_distM, 
                    N = n.ho, M = M, phi = phi, deltasq = deltasq) 
    D_ho <- AD_ho[, M + 1]
    
    C.p <- cbind(X.ho, sparseMatrix(ind_C.ho_X, ind_C.ho_Y,
                          x = c(t(AD_ho[, -(M + 1)]))))
    if(dim(C.p)[2] < length(theta_hat)){theta_hat <- theta_hat[1:dim(C.p)[2]]}
    Y.ho.m <- C.p %*% theta_hat
    
    Y.ho.p <- Y.ho.m + rnorm(n.ho) * sqrt(sigmasq * D_ho)
    
    k.fold.scores[j, "crps"] <- k.fold.scores[j, "crps"] + 
       crps(Y.ho, Y.ho.p@x, sigmasq * D_ho)

    k.fold.scores[j, "rmspe"] <- k.fold.scores[j, "rmspe"] + 
      rmspe(Y.ho, Y.ho.m@x)
    
    cat("\n", i, "th folder,", j, "th phi and deltasq, use time ", 
        (proc.time() - t)[3])
  }
}

save(file = "./data/data/CVscores2.RData", 
          list = c("k.fold.scores"))

