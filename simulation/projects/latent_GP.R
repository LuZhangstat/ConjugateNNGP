
setwd("./simulation")
rm(list = ls())
load("./data/sorted_x_order/nngp_10.RData")


library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
data <- list(N = N, P = P, Y = Y, X = X, dist = c(D_lit),
             as = as, bs = bs, at = at, bt = bt, ap = ap, bp = bp)


# starting value
myinits <-list(list(beta = c(0, 0), sigmasq = 2, tausq = 1, phi = 10, w = w_fit), 
               list(beta = c(1, 1), sigmasq = 3, tausq = 0.5, phi = 6, w = w_fit),
               list(beta = c(-1, -1), sigmasq = 1, tausq = 2, phi = 10, w = w_fit))

parameters <- c("beta", "sigmasq", "tausq", "phi", "w")

samples_w_f <- stan(
  file = "./src/latent_GP.stan",
  data = data,
  init = myinits,
  pars = parameters,
  iter = 2000,
  chains = 3, 
  thin = 1,
  seed = 1
)

# save data in results
save(samples_w_f, file = "./results/latent_GP")

