## draw pic ##

setwd("/Users/luzhang/Documents/github/ConjugateNNGP") # set to the path of ConjugateNNGP
setwd("./SST_study")
rm(list = ls())
library(coda)
library(spBayes)
library(MBA)
library(fields)
library(classInt)
library(RColorBrewer)
library(sp)
library(rworldmap)
library(RgoogleMaps)
library(rworldxtra)
######
library(geoR)
library(raster)
library(leaflet)

# load results

load("./results/conj_LNNGP_paral.RData") 
load("./results/RMSPEresult.RData")
load("./data/data/SST_test.RData")
load("./data/data/SST_data_als.RData")
load("./data/data/SST_train.RData")
load("./data/buildNN/nngp_10.RData")


myqnt <- c(0.5, 0.025, 0.975)
col.qnt <- function(X, start, niter, myqnt){
  return (c(quantile(X[start:niter], myqnt)))
}

newmap <- getMap(resolution = "high")


### plot with package raster

r <- raster(ncols = 480, nrows = 240, xmn = -140, xmx = 0, ymn = 0, 
            ymx = 60)
projection(r) <- "+proj=utm +zone=48 +datum=WGS84"
r_land <- rasterize(newmap, r, getCover = FALSE)
land_extent <- extent(r_land)
t <- proc.time()
#r_cover <- rasterize(SSTdata[, c("lon", "lat")], r, field=1)
#r_length <- rasterize(SSTdata[, c("lon", "lat")], r, fun=function(x ,...) length(x))
r_train_true <- rasterize(SST_train[, c("lon", "lat")], r, SST_train$sst, fun=mean)
r_test_true <- rasterize(SST_test[, c("lon", "lat")], r, SST_test$sst, fun=mean)
r_train_postmean <- rasterize(SST_train[NN.matrix$ord, c("lon", "lat")], r, 
                              gamma_hat[-(1:P)], fun = mean)
r_test_postmean <- rasterize(SST_test[, c("lon", "lat")], r, test_post_mean[, 1], fun=mean)
#p <- data.frame(SSTdata[, c("lon", "lat")], name = SSTdata$sst)
#coordinates(p) <- ~lon+lat
#r_mean <- rasterize(p, r, 'name', fun=mean)
proc.time() - t


col.pal <- colorRampPalette(brewer.pal(11,'RdBu')[1:11])
surf.brks <- classIntervals(r_train_true@data@values, 50, 'pretty')$brks

width <- 620
height <- 360
pointsize <- 16

png(paste("pic/interpol_train_true.png", sep=""), 
    width = width, height = height, pointsize = pointsize, family="Courier")
plot(r_train_true, breaks = surf.brks, interpolate = F, legend = F,
     col = rev(col.pal(length(surf.brks) - 1)), colNA = "yellow")
plot(r_land, add = TRUE, col = 'gray', lwd = 1, legend = F)
plot(r_train_true, legend.only=TRUE, col = rev(col.pal(length(surf.brks) - 1)),
     legend.width=1, legend.shrink = 1,
     axis.args=list(at=seq(surf.brks[1], surf.brks[length(surf.brks)], 5),
                    labels=seq(surf.brks[1], surf.brks[length(surf.brks)], 5), 
                    cex.axis=0.6),
     legend.args=list(text = '', side=4, font=2, line=2.5, cex=0.8))
dev.off()
# white, missing data; gray: land..

png(paste("pic/interpol_test_true.png", sep=""), 
    width = width, height = height, pointsize = pointsize, family="Courier")
plot(r_test_true,  breaks = surf.brks, interpolate = T, legend = F,
     col = rev(col.pal(length(surf.brks) - 1)), colNA = "yellow")
plot(r_land, add = TRUE, col = 'gray', lwd = 1, legend = F)
plot(r_train_true, legend.only=TRUE, col = rev(col.pal(length(surf.brks) - 1)),
     legend.width=1, legend.shrink = 1,
     axis.args=list(at=seq(surf.brks[1], surf.brks[length(surf.brks)], 5),
                    labels=seq(surf.brks[1], surf.brks[length(surf.brks)], 5), 
                    cex.axis=0.6),
     legend.args=list(text = '', side=4, font=2, line=2.5, cex=0.8))
dev.off()

png(paste("pic/interpol_train_post_mean.png", sep=""), 
    width = width, height = height, pointsize = pointsize, family="Courier")
plot(r_train_postmean,  breaks = surf.brks, interpolate = T, legend = F,
     col = rev(col.pal(length(surf.brks) - 1)), colNA = "yellow")
plot(r_land, add = TRUE, col = 'gray', lwd = 1, legend = F)
plot(r_train_true, legend.only=TRUE, col = rev(col.pal(length(surf.brks) - 1)),
     legend.width=1, legend.shrink = 1,
     axis.args=list(at=seq(surf.brks[1], surf.brks[length(surf.brks)], 5),
                    labels=seq(surf.brks[1], surf.brks[length(surf.brks)], 5), 
                    cex.axis=0.6),
     legend.args=list(text = '', side=4, font=2, line=2.5, cex=0.8))
dev.off()

png(paste("pic/interpol_test_post_mean.png", sep=""), 
    width = width, height = height, pointsize = pointsize, family="Courier")
plot(r_test_true,  breaks = surf.brks, interpolate = T, legend = F,
     col = rev(col.pal(length(surf.brks) - 1)), colNA = "yellow")
plot(r_land, add = TRUE, col = 'gray', lwd = 1, legend = F)
plot(r_train_true, legend.only=TRUE, col = rev(col.pal(length(surf.brks) - 1)),
     legend.width=1, legend.shrink = 1,
     axis.args=list(at=seq(surf.brks[1], surf.brks[length(surf.brks)], 5),
                    labels=seq(surf.brks[1], surf.brks[length(surf.brks)], 5), 
                    cex.axis=0.6),
     legend.args=list(text = '', side=4, font=2, line=2.5, cex=0.8))
dev.off()



