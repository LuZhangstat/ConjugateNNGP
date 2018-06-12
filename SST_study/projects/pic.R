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
N = 2544527; P = 3
load("./results/conj_LNNGP_paral.RData")
Conj_pos <- matrix(unlist(Conj_LNNGP_paral[1:100]), nrow = (P + 1 + N))
Conj_pos <- cbind(Conj_pos, 
                  matrix(unlist(Conj_LNNGP_paral[101:200]), nrow = (P + 1 + N)))
Conj_pos <- cbind(Conj_pos, 
                  matrix(unlist(Conj_LNNGP_paral[201:300]), nrow = (P + 1 + N)))
Conj_pos2 <- matrix(unlist(Conj_LNNGP_paral[201:300]), nrow = (P + 1 + N))
Conj_pos1 <- cbind(Conj_pos[1:(P + 1), ], Conj_pos2[1:(P + 1), ])

myqnt <- c(0.5, 0.025, 0.975)
col.qnt <- function(X, start, niter, myqnt){
  return (c(quantile(X[start:niter], myqnt)))
}

t(round(apply(Conj_pos1, 1, col.qnt, 1, 300, myqnt), 2))
#50%  2.5% 97.5%
#[1,]  3.95  3.94  3.95
#[2,] 31.43 31.28 31.59
#[3,]  0.07  0.05  0.09
#[4,] -3.03 -3.08 -2.99

round(rowMeans(Conj_pos1), 2)
# [1]  3.95 31.43  0.07 -3.03

# summary of the std of w
L <- length(Conj_LNNGP_paral)
N <- length(Conj_LNNGP_paral[[1]]$pos.p.w)
w_var <- rep(0, N)
cw <- c()
for(i in 1: N){
  for (j in 1:L){
    cw[j] <- Conj_LNNGP_paral[[j]]$pos.p.w[i]
  }
  w_var[i] <- var(cw)
}
summary(sqrt(w_var))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.06431 0.08466 0.09430 0.09587 0.10593 0.14006 

rm(list = ls())
load("./data/buildNN/nngp_10.RData")
load("./results/RMSPEresult.RData")
load("./data/data/SST_test.RData")
load("./data/data/SST_data_als.RData")
load("./data/data/SST_train.RData")

newmap <- getMap(resolution = "high")
### plot with package raster

r <- raster(ncols = 480, nrows = 240, xmn = -140, xmx = 0, ymn = 0, 
            ymx = 60)
projection(r) <- "+proj=utm +zone=48 +datum=WGS84"
r_land <- rasterize(newmap, r, getCover = FALSE)
land_extent <- extent(r_land)
t <- proc.time()
r_train_true <- rasterize(SST_train[, c("lon", "lat")], r, SST_train$sst, 
                          fun = mean)
r_test_true <- rasterize(SST_test[, c("lon", "lat")], r, SST_test$sst, 
                         fun = mean)
r_train_postmean <- rasterize(SST_train[NN.matrix$ord, c("lon", "lat")], r, 
                              gamma_hat[-(1:P)], fun = mean)
r_test_postmean <- rasterize(SST_test[, c("lon", "lat")], r, 
                             test_post_mean[, 1], fun = mean)
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



