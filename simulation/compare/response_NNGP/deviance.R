
setwd("./simulation")
rm(list = ls())
library("sp")
library("rstan")
library("Matrix")
load("./data/sorted_x_order/nngp_10.RData")
load("./results/response_NNGP")

D_ord <- dist(m.c$coords.ord)
Cov_true <- sigma.sq * exp(- phi * as.matrix(D_ord)) + tau.sq * diag(N)
chol_Cov_true <- chol(Cov_true)
NNGP_ind_i <- c(rep(2:(M), 1:(M-1)), rep(((M+1):N), each = M))
NNGP_ind_j <- c(t(NN.matrix$nearind)[which(t(NN.matrix$nearind) > 0)])

# get posterior samples
post_sim <- extract(samples)

set.seed(1)
n_run <- 500; L <- length(post_sim$sigmasq)
sampleind <- sample.int(L, n_run)

t <- proc.time()
# calculate score and deviance
score <- rep(0, n_run); deviance <- rep(0, n_run)
simind = 0
for (s in sampleind){
  simind = simind + 1
  beta_sim <- post_sim$beta[s,]
  sigmasq_sim <- post_sim$sigmasq[s]
  tausq_sim <- post_sim$tausq[s]
  phi_sim <- post_sim$phi[s]
  kappa_p_1_sim <- tausq_sim / sigmasq_sim + 1;
  
  y_xb_sim <- m.c$y.ord - m.c$X.ord%*%beta_sim
  xb_xbs <- m.c$X.ord %*%(beta_sim - B)
  chol_yxb_sim <- rep(0, N); chol_yxb_sim[1] <- y_xb_sim[1]
  chol_xb_xbs <- rep(0, N); chol_xb_xbs[1] <- xb_xbs[1]
  
  NNGP_L <- rep(0, M * (N - M) + M * (M - 1) / 2)
  D <- rep(sigmasq_sim, N); D[1] <- sigmasq_sim + tausq_sim
  l_record <- 0

  for (i in 2:N){
    # calculate precision matrix for NNGP
    dim <- ifelse(i < (M + 1), (i - 1), M)
    i_nearCorM <- matrix(0, nrow = dim, ncol = dim)
    i_nearCorV <- rep(0, dim)
    # save L as a long vector and use sparse matrix to calculate matrix..
    
    if(dim == 1){
      i_nearCorM[1, 1] <- kappa_p_1_sim 
    }else{
      h <- 0
      for (j in 1:(dim -1)){
        for (k in (j + 1): dim){
          h <- h + 1
          i_nearCorM[j, k] <- exp( - phi_sim * NN.matrix$neardistM[(i - 1), h])
          i_nearCorM[k, j] <- i_nearCorM[j, k]
        }
      }
      diag(i_nearCorM) <- kappa_p_1_sim 
    }
    
    chol_nearCorML <- t(chol(i_nearCorM))
    i_nearCorV <- exp(- phi_sim * NN.matrix$neardist[(i - 1), 1:dim])
    
    v <- forwardsolve(chol_nearCorML, i_nearCorV) # v = chol(CorM)^-1 * CorV
    
    D[i] = D[1] - sigmasq_sim * sum(v^2)
    
    # save CorV %*% inv(CorM)
    NNGP_L[(l_record + 1):(l_record + dim)] <- backsolve(t(chol_nearCorML), v)
    chol_yxb_sim[i] <-y_xb_sim[i] - t(NNGP_L[(l_record+1):(l_record+dim)]) %*% 
      y_xb_sim[NN.matrix$nearind[(i - 1), 1:dim]]
    chol_xb_xbs[i] <- xb_xbs[i] - t(NNGP_L[(l_record+1):(l_record+dim)]) %*% 
      xb_xbs[NN.matrix$nearind[(i - 1), 1:dim]]
    l_record <- l_record + dim
  }
  cat(simind, "th simulation \n")
  score[simind] <- - sum(log(D)) - sum(chol_yxb_sim^2/D)
  I_L <- sparseMatrix(i = 1:N, j = 1:N, x = 1) -
    sparseMatrix(i = NNGP_ind_i, j = NNGP_ind_j, x = NNGP_L, dims = c(N, N))
  

  inv_Cps_Cq <-  (t(I_L / D) %*% I_L)%*% Cov_true
  deviance[simind] <- sum(diag(inv_Cps_Cq)) - 2*sum(log(diag(chol_Cov_true))) +
    sum(log(D)) + sum(chol_xb_xbs^2/D)- N
  
}
proc.time() - t

save(score, deviance, 
     file = "compare/response_NNGP/score_dev.RData")
