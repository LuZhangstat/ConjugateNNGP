setwd(".")
setwd("./simulation")
rm(list = ls())
load("./data/simdata/simdata.RData")
library(spBayes)

starting <- list("beta" = c(-1, 0), "sigma.sq" = 4, "tau.sq" = 2, "phi" = 30)
tuning <- list("phi" = 0.1, "sigma.sq" = 0.1, "tau.sq" = 0.1)
priors <- list("phi.Unif" = c(ap, bp), "sigma.sq.IG" = c(as, bs),
               "tau.sq.IG" = c(at, bt))
cov.model <- "exponential"
n.samples <- 20000
n.report <- 2000
verbose <- TRUE

set.seed(1234)
t <- proc.time()
m.fullGP <- spLM(Y ~ X - 1, coords = coords, starting = starting,
                 tuning = tuning, priors = priors, cov.model = cov.model,
                 n.samples = n.samples, verbose = verbose, n.report = n.report)
proc.time() - t
#user   system  elapsed 
#2446.321   13.308 2499.103 
t <- proc.time()
burn.in <- 0.5 * n.samples + 1
##recover beta and spatial random effects
m.fullGP <- spRecover(m.fullGP, start = burn.in, verbose = FALSE)
proc.time() - t
#user    system   elapsed 
#17228.567    60.306 23147.068

plot(m.fullGP$p.theta.samples, density = FALSE)
m.fullGP.w.summary <- 
  summary(mcmc(t(m.fullGP$p.w.recover.samples)))$quantiles[, c(3, 1, 5)]
plot(w, m.fullGP.w.summary[,1], xlab = "Observed w", ylab = "Fitted w",
     xlim = range(w), ylim = range(m.fullGP.w.summary), 
     main = "Spatial random effects")
arrows(w, m.fullGP.w.summary[, 1], w, m.fullGP.w.summary[, 2], 
       length = 0.02, angle = 90)
lines(range(w), range(w))

save(m.fullGP, file = "./results/fullGP.RData")

