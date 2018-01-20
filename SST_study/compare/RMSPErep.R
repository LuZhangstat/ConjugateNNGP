setwd("") # set to the path of ConjugateNNGP
setwd("./SST_study") 
rm(list = ls())
library("sp")
library("rstan")
library("Matrix")
library("fields")

load("./results/response_NNGP") 
load("./data/data/SST_test.RData")
load("./data/buildNN/nngp_10.RData")

set.seed(1)
# get posterior samples
post_sim <- extract(samples)
n_run <- 5; L <- length(post_sim$sigmasq)
sampleind <- sample.int(L, n_run)

N_test <- dim(SST_test)[1]
coords_test <- cbind(SST_test$projX, SST_test$projY)
X_test <- as.matrix(cbind(1, SST_test[, c("projX", "projY")]))
Y_pred <- Matrix(0, nrow = N_test, ncol = n_run)
X <- as.matrix(X)


RMSPEfun2 <- function(s, l, NNind_l, NNdist, NNdistv, NNdistM){
    beta_sim <- post_sim$beta[s, ]
    sigmasq_sim <- post_sim$sigmasq[s]
    tausq_sim <- post_sim$tausq[s]
    phi_sim <- post_sim$phi[s]
    y_xb_sim <- Y[NNind_l] - X[NNind_l, ] %*% beta_sim
    
    NNcorM <- exp(- phi_sim * NNdistM ) + (tausq_sim / sigmasq_sim) * diag(M) 
    chol_NNcorM <- t(chol(NNcorM))
    NNcorV <- exp(- phi_sim * NNdistv)
    v <- forwardsolve(chol_NNcorM, NNcorV)
    
    var_l <- sigmasq_sim * (1 - sum(v^2))
    
    mean_l <- X_test[l,] %*% beta_sim + 
      t(backsolve(t(chol_NNcorM), v)) %*% y_xb_sim
    Y_pred <- mean_l + rnorm(1)* sqrt(var_l)
    return(Y_pred)
 }


RMSPEfun <- function(l){
  cat(l, "th new sample: \n")
  NNind_l <- sort(order(
    rdist(coords_test[c(1, l),], coords[c(1:N), ])[-1, ])[1:M])
  NNdist <- dist(rbind(coords_test[l, ], coords[NNind_l, ]))
  NNdistv <- c(NNdist)[1:M]
  NNdistM <- as.matrix(NNdist)[-1, -1]

  Y_pred <- unlist(lapply(sampleind, RMSPEfun2, l, NNind_l, NNdist, 
                          NNdistv, NNdistM))
  return(Y_pred)
}

t <- proc.time()
RMSPEresult <- sapply(1:N_test, RMSPEfun)
proc.time() - t


MSPE <- mean((SST_test[1:N_test, 1] - colMeans(RMSPEresult[ , 1:N_test]))^2)
RMSPE <- sqrt(MSPE); RMSPE

save(Y_pred, file = "compare/Y_pred_MSPE_M10_resp.RData")

