setwd("/Users/luzhang/Documents/github/ConjugateNNGP") # set to the path of ConjugateNNGP
setwd("SST_study")
rm(list = ls())
load("./data/buildNN/nngp_10.RData")
source("./projects/functions.R")
library(foreach)
library(doParallel)
library(parallel)

## set estimates for phi and deltasq from cross-validation "spConjNNGP" ##
phi <- 7; deltasq <- 0.001

a = as; b = bs;
ind_x <-c(c(rep(2:M, times = 1:(M - 1)), rep(((M + 1) : N), each = M)), 1:N)
ind_y <- c(c(t(NN.matrix$NN_ind))[which(c(t(NN.matrix$NN_ind)) > 0)], 1:N)
X.ord <- X[NN.matrix$ord, ]
Y.ord <- Y[NN.matrix$ord]

## obtain A and D using C and N(i)
AD <- getADstan(neardist = NN.matrix$NN_dist,  neardistM = NN.matrix$NN_distM, 
                N = N, M = M, phi = phi) 
D <- AD[, M + 1]


## generate sparse matrix X** and Y** 

ind_x_X <- rep(1:N, P); ind_y_X <- rep(1:P, each = N)
ind_x_X_up <- c(ind_x_X, 1:N)
ind_y_X_up <- c(ind_y_X, (P + 1):(N + P))

X_star_up <- 
  sparseMatrix(ind_x_X_up, ind_y_X_up, x = c(1 / sqrt(deltasq) * c(X.ord),
                                             rep(1 / sqrt(deltasq), N)))
X_star_down1 <- 
  sparseMatrix(i = ind_x, j = (ind_y + P), 
               x = c( - c(t(AD[, 1:M]))[
                 which(!is.na(t(AD[, 1:M])))], rep(1, N)) )

#Dsqrtinv <- sparseMatrix(i = 1:N, j = 1:N, x = 1 / sqrt(D))
X_star_down <- X_star_down1 / sqrt(D)
X_star <- rbind(X_star_up, X_star_down)
Y_star <- c(Y.ord / sqrt(deltasq), rep(0, N))

## get theta = (X**^T X**)^-1 X**^T y** by conjugate gradient ##
t <- proc.time()
X_Tstar_Y_star <- t(X_star) %*% Y_star
XTX_star <- t(X_star) %*% X_star
gamma_hat <- cgsparse(XTX_star[1:(N + P), 1:(N + P)],
                      X_Tstar_Y_star[1:(N + P)])
proc.time() - t

## obtain a* and b* by equation for the posterior IG(a*, b*) for sigmasq ##

b_star <- b + 0.5 * sum((Y_star - (X_star %*% gamma_hat))^2) 
a_star <- a + 0.5 * N

## storage for posterior samples sigmasq, beta, and w ##
l <- 300

t_inital <- proc.time()
no_cores <- detectCores() - 2
cl <- makeCluster(no_cores, type = "FORK", outfile = "./results/debug.txt")
clusterSetRNGStream(cl = cl, iseed = 1234)

Conj_LNNGP_paral <- parLapply(cl, 1:l, function(i) {
  
  t <- proc.time()
  
  ## sample sigmasq from IG(a*, b*) ##
  sigmasq <- 1 / rgamma(1, shape = a_star, rate = b_star)
  
  #### sample theta given sigmasq, \hat{\theta} and X** ####
  ##i. sample u ~ N(0, sigmasq * I_{2n+p})$ ##
  
  u <- rnorm(2 * N, 0) * sqrt(sigmasq)
  
  ##ii. get v = (X**T X**)^{-1} X**^Tu$ by conjugate gradient
  
  X_star_u <- t(X_star) %*% u
  v <- cgsparse(XTX_star[1: (N + P), 1:(N + P)], X_star_u[1:(N + P)])
  
  ##iii. get posterior sample theta = \hat{\theta} + v ##
  sample_theta <- gamma_hat + v
  
  cat("\n", i, "th posterior sample:  ")
  cat("sigmasq is ", sigmasq)
  cat(" takes time: ", (proc.time() - t)[3])
  cat(" total time: ", (proc.time() - t_inital)[3], "\n")
  
  return(list(pos.p.sigmasq = sigmasq, pos.p.beta = sample_theta[1:P], 
    pos.p.w = sample_theta[(P + 1):(N + P)]))
  
})

stopCluster(cl)

save(file = "./results/conj_LNNGP_paral.RData",
     list = c("Conj_LNNGP_paral", "gamma_hat"))


