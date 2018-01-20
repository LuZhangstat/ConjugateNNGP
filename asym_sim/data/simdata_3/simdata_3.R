setwd("") # set to the path of ConjugateNNGP
setwd("./asym_sim/data/simdata_3")
rm(list = ls())
library(fields)
library(spBayes)

#---------------------------- Generate simdata --------------------------------#

rmvn <- function(N, mu = 0, V = matrix(1)){
  P <- length(mu)
  if(any(is.na(match(dim(V), P))))
    stop("Dimension problem!")
  D <- chol(V)
  t(matrix(rnorm(N * P), ncol = P) %*% D + rep(mu, rep(N, P)))
}

set.seed(1234)
N <- 900
coords <- cbind(runif(N), runif(N))
X <- as.matrix(cbind(1, rnorm(N)))

B <- as.matrix(c(1, -5))
sigma.sq <- 2
tau.sq <- 0.2
phi <- 30

D <- as.matrix(dist(coords))
R <- exp(-phi*D)
w <- rmvn(1, rep(0, N), sigma.sq*R)
Y <- rnorm(N, X %*% B + w, sqrt(tau.sq))


#----------------------------------- EDA --------------------------------------#
lm.obj <- lm(Y ~ X[, 2])
summary(lm.obj)

library(geoR)
d.max <- max(iDist(coords))
d.max

v.resid <- variog(coords = coords, data = resid(lm.obj), 
                  uvec=(seq(0, 0.4 * d.max, length = 50)))
par(mfrow = c(1, 1))
plot(v.resid, xlab = "Distance (km)")
vario.fit <- variofit(v.resid, cov.model = "exponential")
summary(vario.fit)
plot(v.resid, main = "geoR", cex = 0.1); lines(vario.fit, col = 'red')

#---------------------------------- prior -------------------------------------#
P = 2
at = 2; bt = vario.fit$nugget
as = 2; bs = vario.fit$cov.pars[1]
ap = 3 / d.max ; bp = 3 / (0.01*d.max)
w_fit = rep(0.1, N)
mean_w_zero = rep(0, N)
variofitphi <- 1 / vario.fit$cov.pars[2]

save(list = ls(all.names = TRUE), file = "simdata_3.RData", envir = .GlobalEnv)



