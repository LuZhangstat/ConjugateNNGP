setwd("") # set to the path of ConjugateNNGP
setwd("SST_study")
rm(list = ls())
load("./data/buildNN/nngp_10.RData")

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
data <- list(N = N, M = M, P = P, Y = m.c$y.ord, X = m.c$X.ord, 
             nearind = NN.matrix$nearind, neardist = NN.matrix$neardist, 
             neardistM = NN.matrix$neardistM,
             as = as, bs = bs, at = at, bt = bt, ap = ap, bp = bp)

# starting value
myinits <-list(list(beta = m.c$beta.hat, sigmasq = bs, tausq = bt, 
                    phi = variofitphi), 
               list(beta = c(1, 1, 1), sigmasq = 3, tausq = 0.5, phi = 2),
               list(beta = c(-1, -1, -1), sigmasq = 1, tausq = 0.1, phi = 3))

parameters <- c("beta", "sigmasq", "tausq", "phi")
samples <- stan(
  file = "./src/response_NNGP.stan",
  data = data,
  init = myinits,
  pars = parameters,
  iter = 500, 
  chains = 3, 
  thin = 1,
  seed = 1
)

save(samples, file = "./results/response_NNGP")
