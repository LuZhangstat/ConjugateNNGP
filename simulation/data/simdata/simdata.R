setwd("/Users/luzhang/Documents/github/ConjugateNNGP")
setwd("./simulation/data/simdata")
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
N_total <- 1200
coords_total <- cbind(runif(N_total), runif(N_total))
X_total <- as.matrix(cbind(1, rnorm(N_total)))

B <- as.matrix(c(1, -5))
sigma.sq <- 2
tau.sq <- 0.2
phi <- 30    

D <- as.matrix(dist(coords_total))
R <- exp(-phi*D)
w_total <- rmvn(1, rep(0, N_total), sigma.sq * R)
Y_total <- rnorm(N_total, X_total %*% B + w_total, sqrt(tau.sq))

N = 1000
Y <- Y_total[1:N]
X <- X_total[1:N, ]
w <- w_total[1:N]
coords <- coords_total[1:N, ]
D_lit <- dist(coords)

Y_test <- Y_total[(N + 1) : N_total]
X_test <- X_total[(N + 1) : N_total, ]
w_test <- w_total[(N + 1) : N_total]
coords_test <- coords_total[(N + 1) : N_total, ]


#----------------------------------- EDA --------------------------------------#
lm.obj <- lm(Y ~ X[, 2])
summary(lm.obj)

library(geoR)
d.max <- max(iDist(coords))
d.max

v.resid <- variog(coords = coords, data = resid(lm.obj), 
                  uvec = (seq(0, 0.4*d.max, length = 200)))
par(mfrow=c(1,1))
plot(v.resid, xlab="Distance (km)")
vario.fit <- variofit(v.resid, cov.model="exponential")
summary(vario.fit)
plot(v.resid, main = 'geoR', cex = 0.1); lines(vario.fit, col='red')


#---------------------------------- prior -------------------------------------#
P = 2
at = 2; bt = vario.fit$nugget
as = 2; bs = vario.fit$cov.pars[1]
ap = 3 / d.max ; bp = 3 / (0.01 * d.max)
w_fit = rep(0.1, N)
mean_w_zero = rep(0, N)
variofitphi <- 1 / vario.fit$cov.pars[2]


save(list = ls(all.names = TRUE), file = "simdata.RData", envir = .GlobalEnv)

