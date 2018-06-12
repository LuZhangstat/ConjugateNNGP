setwd("./") # set to the path of ConjugateNNGP
setwd("./simulation")
rm(list = ls())
library("sp")
library("Matrix")
library("fields")
library("RANN")
source("./src/NNMatrix.R")
source("./src/functions.R")
source("./src/util.R")

load("./results/Conj_LNNGP_paral.RData") 
load("./data/simdata/simdata.RData")
load("./data/buildNN/nngp_10.RData")

## set estimates for phi and deltasq from cross-validation "spConjNNGP" ##
phi <- 17.64706; deltasq <- 0.09000

set.seed(1)
# get posterior samples
n_run <- 300; L <- length(Conj_LNNGP_paral)
sampleind <- sample.int(L, n_run)

N_test <- dim(X_test)[1]
Y_pred <- Matrix(0, nrow = N_test, ncol = n_run)
w_pred <- Matrix(0, nrow = N_test, ncol = n_run)
X.ord <- X[NN.matrix$ord, ]
Y.ord <- Y[NN.matrix$ord]
coords.ord <- coords[NN.matrix$ord, ]


t <- proc.time()
## Build nearest neighbor for withheld locations
nn.mod.ho <- nn2(coords.ord, coords_test, k = M)
nn.mod.ho$NN_distM <- sapply(1:N_test, ho_dist, 
                             nn.mod.ho$nn.idx, coords.ord)
## Build matrix Au and Du
ind_Au_X <- rep(1:N_test, each = M)
ind_Au_Y <- c(t(nn.mod.ho$nn.idx))

## obtain A and D using C and N(i)
AD_ho <- getADstan_ho(neardist = nn.mod.ho$nn.dists,  
                      neardistM = nn.mod.ho$NN_distM, 
                      N = N_test, M = M, phi = phi, deltasq = deltasq) 
D_ho <- AD_ho[M + 1, ]
Xu_Au <- cbind(X_test, sparseMatrix(ind_Au_X, ind_Au_Y,
                                x = as.vector(AD_ho[-(M + 1), ])))
L.gamma <- dim(Xu_Au)[2]
proc.time() - t

RMSPEfun <- function(s){
    cat(s, "th sample: \n")
    gamma_sim <- c(Conj_LNNGP_paral[[s]]$pos.p.beta, 
                   Conj_LNNGP_paral[[s]]$pos.p.w)
    sigmasq_sim <- Conj_LNNGP_paral[[s]]$pos.p.sigmasq
    tausq_sim <- sigmasq_sim * deltasq

    Y.ho.m <- Xu_Au %*% gamma_sim[1:L.gamma]
    Y.ho.p <- Y.ho.m + rnorm(N_test) * sqrt(sigmasq_sim * D_ho)
    return(list(Y_pred = Y.ho.p@x))
}

RMSPEresult <- Matrix(unlist(sapply(sampleind, RMSPEfun)), ncol = n_run)


MSPE <- mean((Y_test - rowMeans(RMSPEresult))^2)
RMSPE <- sqrt(MSPE); RMSPE

test_post_mean <- Xu_Au %*% gamma_hat[1:L.gamma]

save(file = "./compare/Conj_NNGP/RMSPEresult.RData",
     list = c("RMSPEresult", "test_post_mean"))





