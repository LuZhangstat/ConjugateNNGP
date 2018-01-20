
setwd("./simulation")
rm(list = ls())
load("data/simdata/simdata.RData")
library(spNNGP)
library(MBA)
library(fields)

##Fit a Conjugate NNGP model and predict for the holdout
sigma.sq.IG <- c(as, bs)

cov.model <- "exponential"

theta.alpha <- as.matrix(expand.grid(3 / seq(0.05, 1, by = 0.05), 
                                     seq(0.1, 2, by = 0.1)))
dim(theta.alpha)

colnames(theta.alpha) <- c("phi", "alpha")
set.seed(1234)
m.c1 <- spConjNNGP(Y~X[, 2], coords=coords, n.neighbors = 10,
                  k.fold = 5, score.rule = "rmspe",
                  n.omp.threads = 2,
                  theta.alpha = theta.alpha, sigma.sq.IG = sigma.sq.IG, 
                  cov.model = cov.model)

theta.alpha <- as.matrix(expand.grid(seq(25, 45, by = 1), seq(0, 0.2,by=0.01)))
dim(theta.alpha)

colnames(theta.alpha) <- c("phi", "alpha")

set.seed(1234)
m.c2 <- spConjNNGP(Y~X[, 2], coords=coords, n.neighbors = 10,
                   k.fold = 5, score.rule = "rmspe",
                   n.omp.threads = 2,
                   theta.alpha = theta.alpha, sigma.sq.IG = sigma.sq.IG, 
                   cov.model = cov.model)
 

m.c2$theta.alpha

width <- 5
height <- 5
pointsize <- 16

pdf(paste("pic/spNNGPCVsim.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")
crps.surf <- mba.surf(m.c2$k.fold.scores[,c("phi","alpha","rmspe")], no.X=100, no.Y=100)$xyz.est
image.plot(crps.surf, xlab=expression(phi), ylab=expression(alpha=tau^2/sigma^2)) 
points(phi, tau.sq/sigma.sq, col="white", pch=0)
points(m.c2$theta.alpha, col="white", pch=23)
legend("topright", legend=c("True", "RMSPE minimum"), col="white", pch=c(0, 23),
       bg = "white", bty="n", text.col = "black")
dev.off()

m.c1$beta.hat; m.c2$beta.hat

m.c1$sigma.sq.hat; m.c2$sigma.sq.hat



