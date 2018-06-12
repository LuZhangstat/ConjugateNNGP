setwd("/Users/luzhang/Documents/github/ConjugateNNGP") # set to the path of ConjugateNNGP
setwd("./simulation")
rm(list = ls())
library(spNNGP)
library(MBA)
library(fields)
load("./data/buildNN/nngp_10.RData")

##Fit a Conjugate NNGP model and predict for the holdout
sigma.sq.IG <- c(as, bs)
cov.model <- "exponential"
theta.alpha <- as.matrix(expand.grid(3 / seq(0.05, 0.5, by = 0.02), 
                                     seq(0.01, 0.2, by = 0.02)))
colnames(theta.alpha) <- c("phi", "alpha")

set.seed(1234)
m.c <- spConjNNGP(Y ~ X-1, coords = coords, 
                  n.neighbors = 10, k.fold = 5, score.rule = "rmspe",
                  n.omp.threads = 2,
                  theta.alpha = theta.alpha, sigma.sq.IG = sigma.sq.IG, 
                  cov.model = cov.model, X.0 = X_test, coords.0 = coords_test)

#> m.c$run.time
#user  system elapsed 
#4.562   0.047   2.329 
#Set phi=17.64706 and alpha=0.09000
round(m.c$beta.hat, 2)
m.c$beta.var
round(m.c$beta.hat[1] - 1.96*sqrt(m.c$beta.var[1, 1]), 2); 
round(m.c$beta.hat[1] + 1.96*sqrt(m.c$beta.var[1, 1]), 2)
round(m.c$beta.hat[2] - 1.96*sqrt(m.c$beta.var[2, 2]), 2); 
round(m.c$beta.hat[2] + 1.96*sqrt(m.c$beta.var[2, 2]), 2)
m.c$theta.alpha

library(invgamma)
a <- m.c$sigma.sq.hat^2 / m.c$sigma.sq.var + 2
b <- m.c$sigma.sq.hat * (a - 1)
round(m.c$sigma.sq.hat, 2)
round(qinvgamma(0.975, a, b), 2)
round(qinvgamma(0.025, a, b), 2)
m.c$k.fold.scores

#RMSPE
MSPE <- mean((Y_test - m.c$y.0.hat)^2)
RMSPE <- sqrt(MSPE); round(RMSPE, 2)

save(file = "./results/Conj_RNNGP.RData",
     list = c("m.c"))

#------------------------- RMSPE interpolated map -----------------------------#
width <- 360
height <- 360
pointsize <- 16

png(paste("pic/spNNGPCV.png", sep = ""), width = width, 
    height = height, pointsize = pointsize, family = "Courier")
ind <- which((m.c$k.fold.scores[, c("alpha")] < 0.1) & 
              m.c$k.fold.scores[, c("phi")])
cvdata <- cbind(m.c$k.fold.scores[ind, c("phi")], 
                m.c$k.fold.scores[ind, c("alpha")], 
                m.c$k.fold.scores[ind, "rmspe"])
cvdata <- m.c$k.fold.scores[, c("phi", "alpha", "rmspe")]
rmspe.surf <- mba.surf(cvdata, no.X = 100, no.Y = 100)$xyz.est
image.plot(rmspe.surf, xlab = expression(phi), 
           ylab = expression(alpha = tau^2 / sigma^2)) 
points(m.c$theta.alpha, col = "white", pch = 2)
legend("topright", legend = c("RMSPE minimum"), col = "white", 
       pch = c(2), bty = "n")

dev.off()


