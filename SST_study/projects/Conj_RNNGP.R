setwd("/Users/luzhang/Documents/github/ConjugateNNGP") # set to the path of ConjugateNNGP
setwd("./SST_study")
rm(list = ls())
library(spNNGP)
library(MBA)
library(fields)
load("./data/buildNN/nngp_10.RData")
load("./data/data/SST_test.RData")

colnames(X) <- c("intercept", "projX", "projY")

sigma.sq.IG <- c(as, bs)
cov.model <- "exponential"
theta.alpha <- as.matrix(expand.grid(seq(3, 21, by = 2), 
                                     seq(0.001, 0.02, by = 0.002)))
set.seed(123)
colnames(theta.alpha) <- c("phi", "alpha")
n.0 <- nrow(SST_test)
m.c <- spConjNNGP(Y ~ X[, "projX"] + X[, "projY"],
                   coords = cbind(X[, "projX"], X[, "projY"]), 
                   n.neighbors = 10, k.fold = 5, score.rule = "rmspe",
                   n.omp.threads = 2, theta.alpha = theta.alpha, 
                   sigma.sq.IG = sigma.sq.IG, cov.model = cov.model, 
                   X.0 = as.matrix(
                     cbind(rep(1, n.0), SST_test[, c("projX", "projY")])), 
                   coords.0 = as.matrix(SST_test[, c("projX", "projY")]))

save(file = "./results/conj_RNNGP.RData",
     list = c("m.c"))
m.c$run.time
round(m.c$beta.hat, 2)
m.c$beta.var
round(m.c$beta.hat[1] - 1.96*sqrt(m.c$beta.var[1, 1]), 2); 
round(m.c$beta.hat[1] + 1.96*sqrt(m.c$beta.var[1, 1]), 2)
round(m.c$beta.hat[2] - 1.96*sqrt(m.c$beta.var[2, 2]), 2); 
round(m.c$beta.hat[2] + 1.96*sqrt(m.c$beta.var[2, 2]), 2)
round(m.c$beta.hat[3] - 1.96*sqrt(m.c$beta.var[3, 3]), 2); 
round(m.c$beta.hat[3] + 1.96*sqrt(m.c$beta.var[3, 3]), 2)
m.c$theta.alpha
round(m.c$sigma.sq.hat, 2)
a <- m.c$sigma.sq.hat^2 / m.c$sigma.sq.var + 2
b <- m.c$sigma.sq.hat * (a - 1)
library(invgamma)
round(qinvgamma(0.975, a, b), 2)
round(qinvgamma(0.025, a, b), 2)


#RMSPE
MSPE <- mean((SST_test$sst - m.c$y.0.hat)^2)
RMSPE <- sqrt(MSPE); round(RMSPE, 2)

