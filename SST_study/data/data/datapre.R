setwd("/Users/luzhang/Documents/github/ConjugateNNGP") # set to the path of ConjugateNNGP
setwd("./SST_study")
rm(list = ls())

## ---------------------------- read in the data ---------------------------- ##

load("./data/data/SST_data_als.RData")

library(mapproj)
knots.sinusoidal <- mapproject(SSTdata$lon, SSTdata$lat, 
                               projection = "sinusoidal")

radius.of.earth = 6.371            ## 6.371 * 1000 kilometers 
knots.sinusoidal = radius.of.earth * (cbind(knots.sinusoidal$x, 
                                            knots.sinusoidal$y))
SSTdata$projX <- knots.sinusoidal[, 1]
SSTdata$projY <- knots.sinusoidal[, 2]
dim(SSTdata)

hist(SSTdata$sst)

#------------------------- get train and test data -------------------------#

set.seed(1)
test_ind <- sample.int(dim(SSTdata)[1], floor(0.1 * length(SSTdata$sst)))
SST_train <- SSTdata[-test_ind, ]
SST_test <- SSTdata[test_ind, ]

save(SST_train, file = "./data/data/SST_train.RData")
save(SST_test, file = "./data/data/SST_test.RData")

#----------------------------- draw pic to check -----------------------------#
library(coda)
library(spBayes)
library(MBA)
library(fields)
library(classInt)
library(RColorBrewer)
library(sp)

par(mfrow= c(1, 1))
h <- 12
surf.raw <- mba.surf(SST_train[c("lon", "lat", "sst")], no.X = 300, 
                     no.Y = 300, exten = F, sp = TRUE, h = h)$xyz.est

surf.brks <- classIntervals(surf.raw$z, 50, 'pretty')$brks
col.pal <- colorRampPalette(brewer.pal(11, 'RdBu')[1:11])
xlim <- range(SST_train$lon)
zlim <- range(surf.raw[["z"]][which(!is.na(surf.raw[["z"]]))])

width <- 10
height <- 8
pointsize <- 16


# plot data on worldmap (for fun)#
# library(rworldmap)
# col.pal <- colorRampPalette(brewer.pal(11,'RdBu')[1:11])
# newmap <- getMap(resolution = "low")
# plot(newmap, xlim = c(-140, 0), ylim = c(0, 60), asp = 1)
# colors <- rev(col.pal(101))
# zcolor <- colors[(SST_test$sst - min(SST_test$sst)) /
#                    diff(range(SST_test$sst))*100 + 1]
# points(SST_test$lon, SST_test$lat, col = zcolor, cex = .1)



#----------------------------------- EDA --------------------------------------#
rm(list = ls())
load("./data/data/SST_train.RData")
Y <- SST_train$sst
X <- cbind(1, SST_train$projX, SST_train$projY)
lm.obj <- lm(Y ~ X[, c(2, 3)])
summary(lm.obj)
coords <- X[, c(2, 3)]
N <- dim(SST_train)[1]
set.seed(123)
subind <- sample.int(N, round(N*0.01))

# Check the variogram 
library(geoR)
d.max <- sqrt((max(SST_train$projX) - min(SST_train$projX))^2 + 
  (max(SST_train$projY) - min(SST_train$projY))^2)
d.max # around 17,000 KM

v.resid <- variog(coords = coords[subind, ], data = resid(lm.obj)[subind], 
                  uvec = (seq(0, 0.5 * d.max, length = 30))) # 30

par(mfrow=c(1,1))
vario.fit <- variofit(v.resid, cov.model="exponential")
summary(vario.fit)

plot(v.resid, xlab = "Distance (1000km)", cex = 0.1)
lines(vario.fit, col='red') 
abline(h = c(vario.fit$nugget, vario.fit$nugget + vario.fit$cov.pars[1]))
abline(v = 3 * vario.fit$cov.pars[2])


#---------------------------------- prior -------------------------------------#
P = 3
at = 2; bt = 0.001 * vario.fit$cov.pars[1]
as = 2; bs = vario.fit$cov.pars[1]
ap = 3 / d.max ; bp = 3 / (0.01 * d.max)
variofitphi <- 1 / vario.fit$cov.pars[2]

save(list = ls(all.names = TRUE), file = "./data/data/realdata.RData", 
     envir = .GlobalEnv)


#---------------------------fit Bayesian linear model--------------------------#
rm(list = ls())
library(spBayes)
load("./data/data/realdata.RData")
n <- nrow(X)
p <- 3

set.seed(1234)
t <- proc.time()
m.1 <- 
  bayesLMConjugate(Y~X[, c(2, 3)], n.samples = 300, 
                   beta.prior.mean = rep(0, times = p),
                   beta.prior.precision = matrix(0, nrow=p, ncol=p),
                   prior.shape = 2, prior.rate = 1)
round(summary(m.1$p.beta.tauSq.samples)$statistics, 2)
round(summary(m.1$p.beta.tauSq.samples)$quantiles, 2)
proc.time() - t

## posterior predictive process 
load("./data/data/SST_test.RData")
n.samples <- 300
coords.test <- cbind(SST_test$projX, SST_test$projY)
X.test <- cbind(1, coords.test)
m.1.pred <- spPredict(m.1, pred.covars = X.test, pred.coords = coords.test,
                      start = 1)

y.hat <- apply(m.1.pred$p.y.predictive.samples, 1, mean)
y.test <- SST_test$sst

RMSPE <- sqrt(sum((y.test - y.hat)^2) / n); round(RMSPE, 2)

print(object.size(x = lapply(ls(), get)), units = "Mb") 
