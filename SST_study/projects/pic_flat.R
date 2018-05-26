# ---- draw pic ---- #
setwd("H:/research/proj1/SST_lu4")
rm(list=ls())

library(rstan)
#library(rstanarm)
load("./data/sorted_x_order/nngp_10.RData")
load("./data/data/SST_train.RData")
load("./data/data/SST_test.RData")

coord2 <- SST_train[m.c$ord, c("lon", "lat")] 


myqnt <- c(0.5, 0.025, 0.975)
col.qnt <- function(X, start, niter, myqnt){
   return (c(quantile(X[start:niter], myqnt)))
 }


load("./results/simpostw_sortx10_conj.RData")
w_conj_stan.qnt <-apply(pos.p.w, 1, col.qnt, 
                   1, 750, myqnt)

library("ggplot2")
library("gridExtra")
#library("TeachingDemos")  #subplot
#library("cowplot")
#library("grid")

#--------------------- fitted w compare ---------------------#
library(coda)
library(spBayes)
library(MBA)
library(fields)
library(classInt)
library(RColorBrewer)
library(sp)

h <- 12

library(rworldmap)
library(RgoogleMaps)
library(rworldxtra)
col.pal <- colorRampPalette(brewer.pal(11,'RdBu')[1:11])
newmap <- getMap(resolution = "high")
#plot(newmap, xlim = c(120, 140), ylim = c(20, 60), asp = 1)

width <- 10
height <- 8
pointsize <- 16
pdf(paste("./pic/traindatamap.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")
xlim <- range(coord2$lon)

plot(newmap, xlim = c(-140, -115), ylim = c(30, 60), asp = 0.58,
     ylab="latitude", xlab="longitude")
axis(2, las=1)
axis(1)
#points(SST_test$lon, SST_test$lat, col = "red", cex = .1)
colors <- rev(col.pal(101)) 
zcolor <- colors[(m.c$y.ord - min(m.c$y.ord)) /
                   diff(range(m.c$y.ord))*100 + 1]
points(coord2$lon, coord2$lat, col = zcolor, cex = .1)
image.plot(legend.only=T, zlim=range(m.c$y.ord), col=colors, 
           legend.cex = 6, legend.mar = 4)
dev.off()


par(mfrow= c(1, 1))
h <- 12
surf.test <- mba.surf(SST_test[c("lon", "lat", "sst")], no.X=300, 
                     no.Y=300, exten=F, sp=TRUE, h=h)$xyz.est

surf.brks <- classIntervals(surf.test$z, 50, 'pretty')$brks
col.pal <- colorRampPalette(brewer.pal(11,'RdBu')[1:11])
#xlim <- c(-10, 50)
xlim <- range(SST_test$lon)
zlim <- range(surf.test[["z"]][which(!is.na(surf.test[["z"]]))])

pdf(paste("./pic/testrawmap.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")

i <- as.image.SpatialGridDataFrame(surf.test)
plot(newmap, xlim = c(-140, -115), ylim = c(30, 60), asp = 0.58,
     ylab="latitude", xlab="longitude")
axis(2, las=1)
axis(1)
image.plot(i, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlim,
           border = "grey50", lwd = 3, legend.mar = 4)
dev.off()


plot(SST_train[c("lon", "lat")], cex = 0.1)




width <- 10
height <- 8
pointsize <- 16
pdf(paste("./pic/conjw10.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")
xlim <- range(coord2$lon)

plot(newmap, xlim = c(-140, -115), ylim = c(30, 60), asp = 0.58,
     ylab="latitude", xlab="longitude")
axis(2, las=1)
axis(1)
#points(SST_test$lon, SST_test$lat, col = "red", cex = .1)
colors <- rev(col.pal(101)) 
zcolor <- colors[(w_conj_stan.qnt[1, ] - min(w_conj_stan.qnt[1, ])) /
                   diff(range(w_conj_stan.qnt[1, ]))*100 + 1]
points(coord2$lon, coord2$lat, col = zcolor, cex = .1)
image.plot(legend.only=T, zlim=range(w_conj_stan.qnt[1, ]), col=colors, 
           legend.cex = 6, legend.mar = 4)
dev.off()



set.seed(1234)
load("./results/simpostw_sortx10_pseudo.RData")
w_pos_stan.qnt <-apply(pos.w, 1, col.qnt, 
                   1, 750, myqnt)

width <- 10
height <- 8
pointsize <- 16
pdf(paste("./pic/pseudow10.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")
xlim <- range(coord2$lon)

plot(newmap, xlim = c(-140, -115), ylim = c(30, 60), asp = 0.58,
     ylab="latitude", xlab="longitude")
axis(2, las=1)
axis(1)
#points(SST_test$lon, SST_test$lat, col = "red", cex = .1)
colors <- rev(col.pal(101)) 
zcolor <- colors[(w_pos_stan.qnt[1, ] - min(w_pos_stan.qnt[1, ])) /
                   diff(range(w_pos_stan.qnt[1, ]))*100 + 1]
points(coord2$lon, coord2$lat, col = zcolor, cex = .1)
image.plot(legend.only=T, zlim=range(w_pos_stan.qnt[1, ]), col=colors, 
           legend.cex = 6, legend.mar = 4)
dev.off()

#--------------------- fitted Y compare ---------------------#

## conjugate ##
load("./compare/Y_pred_MSPE10_conj.RData")
load("./data/data/SST_test.RData")

N_test <- dim(SST_test)[1]
coords_test <- SST_test[, c("lon", "lat")]

fit_SST <- RMSPEresult[seq(2, 2*N_test, by = 2),]
Y_conj_fit.qnt <-apply(fit_SST, 1, col.qnt, 
                   1, 500, myqnt)

par(mfrow= c(1, 1))
h <- 12
surf.fit.conj <- mba.surf(cbind(SST_test[c("lon", "lat")], Y_conj_fit.qnt[1, ]), 
no.X=300, no.Y=300, exten=F, sp=TRUE, h=h)$xyz.est

surf.brks <- classIntervals(surf.test$z, 50, 'pretty')$brks
col.pal <- colorRampPalette(brewer.pal(11,'RdBu')[1:11])
xlim <- range(SST_test$lon)
zlim <- range(surf.test[["z"]][which(!is.na(surf.test[["z"]]))])

pdf(paste("./pic/fitYw10conj.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")

j <- as.image.SpatialGridDataFrame(surf.fit.conj)
plot(newmap, xlim = c(-140, -115), ylim = c(30, 60), asp = 0.58,
     ylab="latitude", xlab="longitude")

axis(2, las=1)
axis(1)
image.plot(j, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlim,
           border = "grey50", lwd = 3, legend.mar = 4)
dev.off()

## pseudo sampler ##
load("./compare/Y_pred_MSPE10_pseudo.RData")
N_test <- dim(SST_test)[1]
coords_test <- SST_test[, c("lon", "lat")]

fit_SST <- RMSPEresult[seq(2, 2*N_test, by = 2),]
Y_pseudo_fit.qnt <-apply(fit_SST, 1, col.qnt, 
                   1, 500, myqnt)

par(mfrow= c(1, 1))
h <- 12
surf.fit.pseudo <- mba.surf(cbind(SST_test[c("lon", "lat")], Y_pseudo_fit.qnt[1, ]), 
no.X=300, no.Y=300, exten=F, sp=TRUE, h=h)$xyz.est

surf.brks <- classIntervals(surf.test$z, 50, 'pretty')$brks
col.pal <- colorRampPalette(brewer.pal(11,'RdBu')[1:11])
xlim <- range(SST_test$lon)
zlim <- range(surf.test[["z"]][which(!is.na(surf.test[["z"]]))])

pdf(paste("./pic/fitYw10pseudo.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")

k <- as.image.SpatialGridDataFrame(surf.fit.pseudo)
plot(newmap, xlim = c(-140, -115), ylim = c(30, 60), asp = 0.58,
     ylab="latitude", xlab="longitude")

axis(2, las=1)
axis(1)
image.plot(j, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlim,
           border = "grey50", lwd = 3, legend.mar = 4)
dev.off()


#----------------- summary of conj and pseudo -------------------#
myqnt <- c(0.5, 0.025, 0.975)
quantile(pos.p.beta[1, ], myqnt); mean(pos.p.beta[1, ])
quantile(pos.p.beta[2, ], myqnt); mean(pos.p.beta[2, ])
quantile(pos.p.beta[3, ], myqnt); mean(pos.p.beta[3, ])
quantile(pos.p.sigmasq, myqnt); mean(pos.p.sigmasq)
quantile(pos.p.sigmasq, myqnt)*0.001; mean(pos.p.sigmasq)*0.001






