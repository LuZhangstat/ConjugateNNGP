# ---- draw pic ---- #
setwd("./simulation")
rm(list = ls())
library(spNNGP)
library(spBayes)

load("./data/buildNN/nngp_10.RData")

# load output from all models
load("./results/Conj_LNNGP_paral.RData")
Conj_pos <- matrix(unlist(Conj_LNNGP_paral), nrow = (P + 1 + N))
rm(list = "Conj_LNNGP_paral")
load("./results/fullGP.RData")
load("./results/NNGP_latent.RData")

# inference of process parameters
n.samples = 20000
burn.in <- 0.5 * n.samples + 1
round(summary(m.fullGP$p.beta.recover.samples)$quantiles[,c(3,1,5)],2)
round(summary(window(m.NNGP.s$p.beta.samples, start=burn.in))$quantiles[,c(3,1,5)],2)
round(summary(window(m.fullGP$p.theta.samples, start=burn.in))$quantiles[,c(3,1,5)],2)
round(summary(window(m.NNGP.s$p.theta.samples, start=burn.in))$quantiles[,c(3,1,5)],2)

myqnt <- c(0.5, 0.025, 0.975)
col.qnt <- function(X, start, niter, myqnt){
  return (c(quantile(X[start:niter], myqnt)))
}
t(round(apply(Conj_pos[1:(P + 1), ], 1, col.qnt, 1, 300, myqnt), 2))

# running time:
m.fullGP$run.time
m.NNGP.s$run.time

w_CLNNGP.qnt <- apply(Conj_pos[-(1:(P + 1)), ], 1, col.qnt, 1, 300, myqnt)
w_LNNGP.qnt <- apply(m.NNGP.s$p.w.samples, 1, col.qnt, 
                   burn.in, n.samples, myqnt)
w_LNNGP.m <- rowMeans(m.NNGP.s$p.w.samples[, burn.in:n.samples]) 
w_fullGP.qnt <- 
  t(summary(mcmc(t(m.fullGP$p.w.recover.samples)))$quantiles[, c(3, 1, 5)])
w_fullGP.m <- rowMeans(m.fullGP$p.w.recover.samples) 

plot(w[NN.matrix$ord],  w_CLNNGP.qnt[1, ])
plot(w,  w_LNNGP.qnt[1, ])
plot(w,  w_fullGP.m)
library("ggplot2")
library("gridExtra")
library("grid")

#--------------------- fitted w compare ---------------------#
library(coda)
library(spBayes)
library(MBA)
library(fields)
library(classInt)
library(RColorBrewer)
library(sp)

h <- 12
surf.raw <- mba.surf(cbind(coords, w), no.X = 300, no.Y = 300, 
                     exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.LNNGP <- mba.surf(cbind(coords, w_LNNGP.m), no.X=300, no.Y=300, 
                       exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.CLNNGP <- mba.surf(cbind(coords[NN.matrix$ord, ], gamma_hat[-(1:P)]), 
                        no.X=300, no.Y=300, exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.fullGP <- mba.surf(cbind(coords, w_fullGP.m), 
                        no.X=300, no.Y=300, exten = TRUE, sp = TRUE, h = h)$xyz.est

surf.brks <- classIntervals(surf.raw$z, 500, 'pretty')$brks
col.pal <- colorRampPalette(brewer.pal(11,'RdBu')[11:1])
xlim <- c(0, 1.13)
zlim <- range(c(surf.raw[["z"]], surf.LNNGP[["z"]], surf.CLNNGP[["z"]], 
                surf.fullGP[["z"]]))

# size for the mapping of w               
width <- 360
height <- 360
pointsize <- 16

png(paste("./pic/map-w-true.png", sep = ""), 
    width = width, height = height, pointsize = pointsize, family = "Courier")
par(mfrow = c(1, 1))
##Obs
i <- as.image.SpatialGridDataFrame(surf.raw)
plot(coords, typ="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="y", xlab="x") 
     #main = "true")
axis(2, las=1)
axis(1)
image.plot(i, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlim)
dev.off()

png(paste("./pic/map-w-LNNGP.png", sep=""), 
    width = width, height = height, pointsize = pointsize, family = "Courier")

##LNNGP
j <- as.image.SpatialGridDataFrame(surf.LNNGP)
#pdf(paste("w-obs.pdf", sep=""), width=width, height=height, pointsize=pointsize, family="Courier")
plot(coords, typ = "n", cex = 0.5, xlim = xlim, axes = FALSE, ylab = "y", 
     xlab = "x", main = "")
axis(2, las = 1)
axis(1)
image.plot(j, add = TRUE, col = rev(col.pal(length(surf.brks)-1)), zlim = zlim)
dev.off()

png(paste("./pic/map-w-CLNNGP.png", sep=""),
    width=width, height=height, pointsize=pointsize, family="Courier")
##CLNNGP
k <- as.image.SpatialGridDataFrame(surf.CLNNGP)
plot(coords, typ = "n", cex = 0.5, xlim = xlim, axes = FALSE, ylab="y", 
     xlab="x"#, main = "fullGP"
     )
axis(2, las=1)
axis(1)
image.plot(k, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlim)
dev.off()

png(paste("./pic/map-w-fullGP.png", sep=""),
    width=width, height=height, pointsize=pointsize, family="Courier")
##fullNNGP
l <- as.image.SpatialGridDataFrame(surf.fullGP)
plot(coords, typ = "n", cex = 0.5, xlim = xlim, axes = FALSE, ylab="y", 
     xlab="x")
axis(2, las=1)
axis(1)
image.plot(l, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlim)
dev.off()


#-------------------compare the width of 95% CI------------------#

h <- 12
surf.LNNGPCI <- mba.surf(cbind(coords, (w_LNNGP.qnt[3, ] - w_LNNGP.qnt[2, ])), 
                        no.X = 300, no.Y=300, exten=TRUE, sp=TRUE, h=h)$xyz.est
surf.CLNNGPCI <- mba.surf(
  cbind(coords[NN.matrix$ord, ], (w_CLNNGP.qnt[3, ] - w_CLNNGP.qnt[2, ])),
  no.X = 300, no.Y = 300, exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.fullGPCI <- mba.surf(
  cbind(coords, (w_fullGP.qnt[3, ] - w_fullGP.qnt[2, ])),
  no.X = 300, no.Y = 300, exten = TRUE, sp = TRUE, h = h)$xyz.est

surf.brks <- classIntervals(surf.CLNNGPCI$z, 500, 'pretty')$brks
col.pal <- colorRampPalette(brewer.pal(9,'YlGn')[9:1])
xlim <- c(0,1.13)
zlim <- range(c(surf.LNNGPCI[["z"]], surf.CLNNGPCI[["z"]], 
                surf.fullGPCI[["z"]]))

width <- 360
height <- 360
pointsize <- 16

png(paste("./pic/map-w-LNNGPCI.png", sep = ""), 
    width = width, height = height, pointsize = pointsize, family = "Courier")
par(mfrow = c(1, 1))
##nngp
j <- as.image.SpatialGridDataFrame(surf.LNNGPCI)
plot(coords, typ = "n", cex = 0.5, xlim = xlim, axes = FALSE, 
     ylab = "y", xlab = "x")
axis(2, las = 1)
axis(1)
image.plot(j, add = TRUE, 
           col = rev(col.pal(length(surf.brks) - 1)), zlim = zlim)
dev.off()

png(paste("./pic/map-w-CLNNGPCI.png", sep = ""), 
    width = width, height = height, pointsize = pointsize, family = "Courier")
par(mfrow = c(1, 1))
##nngp
k <- as.image.SpatialGridDataFrame(surf.CLNNGPCI)
plot(coords, typ = "n", cex = 0.5, xlim = xlim, axes = FALSE, 
     ylab = "y", xlab = "x")
axis(2, las = 1)
axis(1)
image.plot(k, add = TRUE, 
           col = rev(col.pal(length(surf.brks) - 1)), zlim = zlim)
dev.off()

png(paste("./pic/map-w-fullGPCI.png", sep = ""), 
    width = width, height = height, pointsize = pointsize, family = "Courier")
par(mfrow = c(1, 1))
##nngp
l <- as.image.SpatialGridDataFrame(surf.fullGPCI)
plot(coords, typ = "n", cex = 0.5, xlim = xlim, axes = FALSE, 
     ylab = "y", xlab = "x")
axis(2, las = 1)
axis(1)
image.plot(l, add = TRUE, 
           col = rev(col.pal(length(surf.brks) - 1)), zlim = zlim)
dev.off()

#-------------------QQ plots------------------#
w_CLNNGP_data <- data.frame(w_true = w[NN.matrix$ord], 
                            w_median = w_CLNNGP.qnt[1, ], 
                            w_2.5 = w_CLNNGP.qnt[2, ], 
                            w_97.5 = w_CLNNGP.qnt[3, ])
w_LNNGP_data <- data.frame(w_true = w, w_median = w_LNNGP.qnt[1, ], 
                           w_2.5 = w_LNNGP.qnt[2, ], 
                           w_97.5 = w_LNNGP.qnt[3, ])
w_fullGP_data <- data.frame(w_true = w, w_median = w_fullGP.qnt[1, ], 
                           w_2.5 = w_fullGP.qnt[2, ], 
                           w_97.5 = w_fullGP.qnt[3, ])

require(ggplot2)
png(paste("./pic/w-CLNNGPCI.png", sep = ""), 
    width = width, height = height, pointsize = pointsize, family = "Courier")
ggplot(w_CLNNGP_data, aes(x = w_true, y = w_median)) +
  geom_point(size = 0.2) + 
  geom_errorbar(aes(ymax = w_97.5, ymin = w_2.5), size = 0.1) +
  geom_abline(intercept = 0, slope = 1)
dev.off()
sum((w_CLNNGP_data$w_true > w_CLNNGP_data$w_2.5) &
      (w_CLNNGP_data$w_true < w_CLNNGP_data$w_97.5))
# 955

png(paste("./pic/w-LNNGPCI.png", sep = ""), 
    width = width, height = height, pointsize = pointsize, family = "Courier")
ggplot(w_LNNGP_data, aes(x = w_true, y = w_median)) +
  geom_point(size = 0.2) + 
  geom_errorbar(aes(ymax = w_97.5, ymin = w_2.5), size = 0.1) +
  geom_abline(intercept = 0, slope = 1)
dev.off()
sum((w_LNNGP_data$w_true > w_LNNGP_data$w_2.5) &
      (w_LNNGP_data$w_true < w_LNNGP_data$w_97.5))
#951

png(paste("./pic/w-fullGPCI.png", sep = ""), 
    width = width, height = height, pointsize = pointsize, family = "Courier")
ggplot(w_fullGP_data, aes(x = w_true, y = w_median)) +
  geom_point(size = 0.2) + 
  geom_errorbar(aes(ymax = w_97.5, ymin = w_2.5), size = 0.1) +
  geom_abline(intercept = 0, slope = 1)
dev.off()
sum((w_fullGP_data$w_true > w_fullGP_data$w_2.5) &
      (w_fullGP_data$w_true < w_fullGP_data$w_97.5))
#946


qqplot(m.NNGP.s$p.w.samples.ord[200, burn.in:n.samples], Conj_pos[(P + 1)+200, ])
abline(a = 0, b = 1)
