setwd("/Users/luzhang/Documents/github/ConjugateNNGP") # set to the path of ConjugateNNGP
setwd("./SST_study")
rm(list = ls())
library(spNNGP)
library(MBA)
library(fields)
load("./data/buildNN/nngp_10.RData")

colnames(X) <- c("intercept", "projX", "projY")

##Fit a Conjugate NNGP model and predict for the holdout
sigma.sq.IG <- c(as, bs)
cov.model <- "exponential"
theta.alpha <- as.matrix(expand.grid(3 / seq(0.01, 0.4, by = 0.02), 
                                     seq(0.01, 0.2, by = 0.02)))
colnames(theta.alpha) <- c("phi", "alpha")

set.seed(1234)
m.c1 <- spConjNNGP(Y ~ X[, "projX"] + X[, "projY"],
                   coords = cbind(X[, "projX"], X[, "projY"]), 
                   n.neighbors = 10, k.fold = 5, score.rule = "rmspe",
                   n.omp.threads = 4,
                   theta.alpha = theta.alpha, sigma.sq.IG = sigma.sq.IG, 
                   cov.model = cov.model)


theta.alpha <- as.matrix(expand.grid(seq(3, 41, by = 2), 
                                     seq(0.001, 0.02, by = 0.001)))
set.seed(1234)
colnames(theta.alpha) <- c("phi", "alpha")
m.c2 <- spConjNNGP(Y ~ X[, "projX"] + X[, "projY"],
                   coords = cbind(X[, "projX"], X[, "projY"]), 
                   n.neighbors = 10, k.fold = 5, score.rule = "rmspe",
                   n.omp.threads = 6, theta.alpha = theta.alpha, 
                   sigma.sq.IG = sigma.sq.IG, cov.model = cov.model)
m.c2$theta.alpha
save(file = "./results/mc2.RData",
list = c("m.c2"))
load("./results/mc2.RData")

#------------------------- RMSPE interpolated map -----------------------------#
width <- 360
height <- 360
pointsize <- 16

esi <- matrix(m.c1$k.fold.scores[which.min(m.c1$k.fold.scores[, "rmspe"]),
                                    c("phi", "alpha")], nrow = 1, ncol = 2)
colnames(esi) <- c("phi", "alpha")
m.c1$k.fold.scores[order(m.c1$k.fold.scores[, "rmspe"])[1:20], 
                   c("phi", "alpha", "rmspe")]


png(paste("pic/spNNGPCV2.png", sep = ""), width = width, 
    height = height, pointsize = pointsize, family = "Courier")
ind <- which((m.c2$k.fold.scores[, c("alpha")] < 0.1) & 
              m.c2$k.fold.scores[, c("phi")])
cvdata <- cbind(m.c2$k.fold.scores[ind, c("phi")], 
                m.c2$k.fold.scores[ind, c("alpha")], 
                m.c2$k.fold.scores[ind, "rmspe"])
cvdata <- m.c2$k.fold.scores[, c("phi", "alpha", "rmspe")]
rmspe.surf <- mba.surf(cvdata, no.X = 100, no.Y = 100)$xyz.est
image.plot(rmspe.surf, xlab = expression(phi), 
           ylab = expression(alpha = tau^2 / sigma^2)) 
points(m.c2$theta.alpha, col = "white", pch = 2)
legend("topright", legend = c("RMSPE minimum"), col = "white", 
       pch = c(2), bty = "n")

dev.off()

m.c2$beta.hat; m.c2$sigma.sq.hat

