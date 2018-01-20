setwd("") # set to the path of ConjugateNNGP
rm(list = ls())
load("./data/simdata_2/nngp_5_2.RData")

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
data <- list(N = N, M = M, P = P, Y = m.c$y.ord, X = m.c$X.ord, 
             nearind = NN.matrix$nearind, neardist = NN.matrix$neardist, 
             neardistM = NN.matrix$neardistM,
             as = as, bs = bs, at = at, bt = bt, ap = ap, bp = bp)

# starting value
myinits <-list(list(beta = c(0, 0), sigmasq = 2, tausq = 1, phi = 10), 
               list(beta = c(1, 1), sigmasq = 3, tausq = 0.5, phi = 6),
               list(beta = c(-1, -1), sigmasq = 1, tausq = 2, phi = 10))

parameters <- c("beta", "sigmasq", "tausq", "phi")
samples <- stan(
  file = "./src/nngp.stan",
  data = data,
  init = myinits,
  pars = parameters,
  iter = 2000, 
  chains = 3, 
  thin = 1,
  seed = 1
)

save(samples, file = "./results/sim2/response_nngp_5")

