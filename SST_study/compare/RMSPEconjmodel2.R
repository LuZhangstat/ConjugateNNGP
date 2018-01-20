setwd("") # set to the path of ConjugateNNGP
setwd("./SST_study") 
rm(list = ls())
library("sp")
library("rstan")
library("Matrix")
library("fields")

load("./results/conj_model2.RData") 
load("./data/data/SST_test.RData")
load("./data/buildNN/nngp_10.RData")

## set estimates for phi and deltasq from cross-validation "spConjNNGP" ##
phi_sim = 14; deltasq = 0.001

set.seed(1)
# get posterior samples
n_run <- 500; L <- length(pos.p.sigmasq)
sampleind <- sample.int(L, n_run)

N_test <- dim(SST_test)[1]
coords_test <- cbind(SST_test$projX, SST_test$projY)
X_test <- as.matrix(cbind(1, SST_test[, c("projX", "projY")]))
Y_pred <- Matrix(0, nrow = N_test, ncol = n_run)
w_pred <- Matrix(0, nrow = N_test, ncol = n_run)
X.ord <- as.matrix(m.c$X.ord)


RMSPEfun2 <- function(s, l, NNind_l, NNdist, NNdistv, NNdistM){
    beta_sim <- pos.p.beta[, s]
    sigmasq_sim <- pos.p.sigmasq[s]
    tausq_sim <- sigmasq_sim*deltasq
    w_sim <- pos.p.w[, s]
    y_xb_sim <- m.c$y.ord[NNind_l] - X.ord[NNind_l, ]%*%beta_sim
    
    NNcorM <- exp(- phi_sim * NNdistM ) 
    chol_NNcorM <- t(chol(NNcorM))
    NNcorV <- exp(- phi_sim * NNdistv)
    v <- forwardsolve(chol_NNcorM, NNcorV)
    
    var_l <- sigmasq_sim * (1 - sum(v^2))
    
    mean_l <-  t(backsolve(t(chol_NNcorM), v)) %*% w_sim[NNind_l]
    w_pred <- mean_l + rnorm(1)* sqrt(var_l)

    Y_pred <- X_test[l,] %*% beta_sim + w_pred + rnorm(1)* sqrt(tausq_sim)
    return(list(w_pred = w_pred, Y_pred = Y_pred))
  }


RMSPEfun <- function(l){
	cat(l, "th new sample: \n")
  NNind_l <- sort(order(
    rdist(coords_test[c(1, l),], m.c$coords.ord)[-1, ])[1:M])
  NNdist <- dist(rbind(coords_test[l, ], m.c$coords.ord[NNind_l, ]))
  NNdistv <- c(NNdist)[1:M]
  NNdistM <- as.matrix(NNdist)[-1, -1]

  result2 <- Matrix(unlist(lapply(sampleind, RMSPEfun2, l, NNind_l, NNdist, 
                                  NNdistv, NNdistM)), nrow = 2)
  return(list(w_pred = result2[1, ], Y_pred = result2[2, ]))
}

t <- proc.time()
RMSPEresult <- Matrix(unlist(sapply(1:N_test, RMSPEfun)), byrow = T, 
                      ncol = n_run)
proc.time() - t

MSPE <- mean((SST_test[1:N_test, 1] - 
                rowMeans(RMSPEresult[seq(2, 2 * N_test, by = 2),]))^2)
RMSPE <- sqrt(MSPE); RMSPE
 
save(RMSPEresult, file = "compare/Y_pred_MSPE10_conj_model2.RData")





