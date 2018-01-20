# set the working directory below:
# setwd("H:/research/proj1/asym_sim")
rm(list = ls())
load("./data/simdata_5/nngp_5_5.RData")

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
data <- list(N = N, M = M, P = P, Y = m.c$y.ord, X = m.c$X.ord, 
             nearind = NN.matrix$nearind, neardist = NN.matrix$neardist, 
             neardistM = NN.matrix$neardistM,
             as = as, bs = bs, at = at, bt = bt, ap = ap, bp = bp)


# starting value
myinits <-list(list(beta = c(0, 0), sigmasq = 2, tausq = 1, phi = 10, w = w_fit), 
               list(beta = c(1, 1), sigmasq = 3, tausq = 0.5, phi = 6, w = w_fit),
               list(beta = c(-1, -1), sigmasq = 1, tausq = 2, phi = 10, w = w_fit))


parameters <- c("beta", "sigmasq", "tausq", "phi", "w")

samples_w_b1 <- stan(

  file = "./src/nngp_w_b1.stan",
  data = data,
  init = myinits,
  pars = parameters,
  iter = 2000, 
  chains = 3, 
  thin = 1,
  seed = 1
)

#stan_trace(samples_w)
save(samples_w_b1, file = "./results/sim5/model2_nngp_5_b1")



