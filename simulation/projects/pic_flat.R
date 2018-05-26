# ---- draw pic ---- #
setwd("./simulation")
rm(list=ls())

library(rstan)
library(rstanarm)
load("./data/sorted_x_order/nngp_10.RData")
load("./results/sorted_x/model2_nngp_10")
sum_m2_post_sim <- summary(samples_w)
sim_w_post <- extract(samples_w)

load("./results/sorted_x/model2_full")
sum_m2_post_sim_f <- summary(samples_w_f)
sim_w_f_post <- extract(samples_w_f)

myqnt <- c(0.5, 0.025, 0.975)
col.qnt <- function(X, start, niter, myqnt){
   return (c(quantile(X[start:niter], myqnt)))
 }

w_stan.qnt <-apply(sim_w_post$w - sim_w_post$beta[, 1], 2, col.qnt, 
                   1, 3000, myqnt)
w_f_stan.qnt <-apply(sim_w_f_post$w - sim_w_f_post$beta[, 1], 2, col.qnt, 
                     1, 3000, myqnt)

load("./results/sorted_x/simpostw_sortx10_pseudo.RData")
w_fast_stan.qnt <-apply(pos.w, 1, col.qnt, 
                   1, 3000, myqnt)

load("./results/sorted_x/simpostw_sortx10_conj.RData")
w_p_fast_stan.qnt <-apply(pos.p.w, 1, col.qnt, 
                        1, 3000, myqnt)

plot(w[m.c$ord],  w_stan.qnt[1, ])
plot(w,  w_f_stan.qnt[1, ])
plot(w[m.c$ord], w_fast_stan.qnt[1, ])
plot(w_stan.qnt[1, ],  w_f_stan.qnt[1, m.c$ord])
plot(w_stan.qnt[1, ],  w_fast_stan.qnt[1, ])
plot(w_f_stan.qnt[1, m.c$ord],  w_fast_stan.qnt[1, ])
library("ggplot2")
library("gridExtra")
library("grid")

w_nngp_plot_data<- data.frame(w_true = w[m.c$ord], w_median = w_stan.qnt[1, ], 
                              w_2.5 = w_stan.qnt[2, ], w_97.5 = w_stan.qnt[3, ])
w_fullGP_plot_data<- data.frame(w_true = w, w_median = w_f_stan.qnt[1, ], 
                              w_2.5 = w_f_stan.qnt[2, ], w_97.5 = w_f_stan.qnt[3, ])

w_fast_plot_data<- data.frame(w_true = w[m.c$ord], w_median = w_fast_stan.qnt[1, ], 
                              w_2.5 = w_fast_stan.qnt[2, ], 
                              w_97.5 = w_fast_stan.qnt[3, ])

w_sim_plot_data<- data.frame(w_true = w_f_stan.qnt[1, m.c$ord], 
                             w_median = w_stan.qnt[1, ], 
                             w_2.5 = w_stan.qnt[2, ], 
                             w_97.5 = w_stan.qnt[3, ])

w_sim2_plot_data<- data.frame(w_true = w_stan.qnt[1, ], 
                              w_median = w_fast_stan.qnt[1, ], 
                              w_2.5 = w_fast_stan.qnt[2, ], 
                              w_97.5 = w_fast_stan.qnt[3, ])


#--------------------- fitted w compare ---------------------#
library(coda)
library(spBayes)
library(MBA)
library(fields)
library(classInt)
library(RColorBrewer)
library(sp)

h <- 12
surf.raw <- mba.surf(cbind(coords, w), no.X=300, no.Y=300, 
                   exten=TRUE, sp=TRUE, h=h)$xyz.est
surf.nngp <- mba.surf(cbind(coords[m.c$ord, ], w_stan.qnt[1, ]), no.X=300, no.Y=300, 
                   exten=TRUE, sp=TRUE, h=h)$xyz.est
surf.fullGP <- mba.surf(cbind(coords, w_f_stan.qnt[1, ]), no.X=300, no.Y=300, 
                      exten=TRUE, sp=TRUE, h=h)$xyz.est
surf.fast <- mba.surf(cbind(coords[m.c$ord, ], w_fast_stan.qnt[1, ]), no.X=300, no.Y=300, 
                        exten=TRUE, sp=TRUE, h=h)$xyz.est
surf.fast.p <- mba.surf(cbind(m.c$coords.ord, w_p_fast_stan.qnt[1, ]), 
                        no.X=300, no.Y=300, 
                      exten=TRUE, sp=TRUE, h=h)$xyz.est

surf.brks <- classIntervals(surf.raw$z,500, 'pretty')$brks
col.pal <- colorRampPalette(brewer.pal(11,'RdBu')[11:1])
xlim <- c(0,1.13)
zlim <- range(c(surf.raw[["z"]], surf.nngp[["z"]],surf.fullGP[["z"]]))

# size for the mapping of w               
width <- 5
height <- 5
pointsize <- 16

pdf(paste("./pic/map-w-true.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")
par(mfrow = c(1, 1))

##Obs
i <- as.image.SpatialGridDataFrame(surf.raw)
plot(coords, typ="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="y", xlab="x") 
     #main = "true")
axis(2, las=1)
axis(1)
image.plot(i, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlim)
dev.off()

pdf(paste("./pic/map-w-nngp.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")

##nngp
j <- as.image.SpatialGridDataFrame(surf.nngp)
#pdf(paste("w-obs.pdf", sep=""), width=width, height=height, pointsize=pointsize, family="Courier")
plot(coords, typ="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="y", 
     xlab="x", main = "")
axis(2, las=1)
axis(1)
image.plot(j, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlim)
dev.off()

pdf(paste("./pic/map-w-fullGP.pdf", sep=""),
    width=width, height=height, pointsize=pointsize, family="Courier")

##fullGP
k <- as.image.SpatialGridDataFrame(surf.fullGP)
#pdf(paste("w-obs.pdf", sep=""), width=width, height=height, pointsize=pointsize, family="Courier")
plot(coords, typ="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="y", 
     xlab="x"#, main = "fullGP"
     )
axis(2, las=1)
axis(1)
image.plot(k, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlim)
dev.off()

pdf(paste("./pic/map-w-pseudo.pdf", sep=""),
    width=width, height=height, pointsize=pointsize, family="Courier")
l <- as.image.SpatialGridDataFrame(surf.fast)
#pdf(paste("w-obs.pdf", sep=""), width=width, height=height, pointsize=pointsize, family="Courier")
plot(coords, typ="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="y", 
     xlab="x"#, main = "Pseudo"
     )
axis(2, las=1)
axis(1)
image.plot(l, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlim)
dev.off()

pdf(paste("./pic/map-w-conjugate.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")
m <- as.image.SpatialGridDataFrame(surf.fast.p)
#pdf(paste("w-obs.pdf", sep=""), width=width, height=height, pointsize=pointsize, family="Courier")
plot(coords, typ="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="y", 
     xlab="x"  #, main = "conjugate"
     )
axis(2, las=1)
axis(1)
image.plot(m, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlim)
dev.off()

library(geoR)
pdf(paste("./pic/sim_em_variog.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")
plot(v.resid, cex = 0.2, ylim = c(0, 3)); lines(vario.fit)
abline(h = c(bt, bs+bt)); abline(v = 3*vario.fit$cov.pars[2])
dev.off()


#-------------------compare the width of 95% CI------------------#

h <- 12
surf.nngpCI <- mba.surf(cbind(coords, (w_stan.qnt[3, ] - w_stan.qnt[1, ])), 
                        no.X=300, no.Y=300, 
                      exten=TRUE, sp=TRUE, h=h)$xyz.est
surf.fullGPCI <- mba.surf(cbind(coords, (w_f_stan.qnt[3, ] - w_f_stan.qnt[1, ])),
                        no.X=300, no.Y=300, 
                        exten=TRUE, sp=TRUE, h=h)$xyz.est
surf.fastCI <- mba.surf(cbind(coords, (w_fast_stan.qnt[3, ] - w_fast_stan.qnt[1, ])),
                      no.X=300, no.Y=300, 
                      exten=TRUE, sp=TRUE, h=h)$xyz.est
surf.fast.pCI <- mba.surf(cbind(coords, (w_p_fast_stan.qnt[3, ] - w_p_fast_stan.qnt[1, ])),
                        no.X=300, no.Y=300, 
                        exten=TRUE, sp=TRUE, h=h)$xyz.est

surf.brks <- classIntervals(surf.nngpCI$z, 500, 'pretty')$brks
col.pal <- colorRampPalette(brewer.pal(11,'RdBu')[11:1])
xlim <- c(0,1.13)
zlim <- range(c(surf.nngpCI[["z"]], surf.fullGPCI[["z"]],
                surf.fastCI[["z"]], surf.fast.pCI[["z"]]))

#width <- 16
#height <- 6
#pointsize <- 16

width <- 10
height <- 10
pointsize <- 16

pdf(paste("stan_test_paper_flat/pic/w-CI-compare3.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")
par(mfrow = c(2, 2))

##nngp
j <- as.image.SpatialGridDataFrame(surf.nngpCI)
#pdf(paste("w-obs.pdf", sep=""), width=width, height=height, pointsize=pointsize, family="Courier")
plot(coords, typ="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="y", 
     xlab="x", main = "NNGPCI")
axis(2, las=1)
axis(1)
image.plot(j, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlim)

##fullGP
k <- as.image.SpatialGridDataFrame(surf.fullGPCI)
#pdf(paste("w-obs.pdf", sep=""), width=width, height=height, pointsize=pointsize, family="Courier")
plot(coords, typ="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="y", 
     xlab="x", main = "fullGPCI")
axis(2, las=1)
axis(1)
image.plot(k, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlim)

l <- as.image.SpatialGridDataFrame(surf.fastCI)
#pdf(paste("w-obs.pdf", sep=""), width=width, height=height, pointsize=pointsize, family="Courier")
plot(coords, typ="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="y", 
     xlab="x", main = "fastCI")
axis(2, las=1)
axis(1)
image.plot(l, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlim)

m <- as.image.SpatialGridDataFrame(surf.fast.pCI)
#pdf(paste("w-obs.pdf", sep=""), width=width, height=height, pointsize=pointsize, family="Courier")
plot(coords, typ="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="y", 
     xlab="x", main = "fast point CI")
axis(2, las=1)
axis(1)
image.plot(m, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlim)


dev.off()


width <- 5
height <- 5
pointsize <- 16

par(mfrow = c(1, 3))

qqplot(sim_w_post$sigma^2, pos.sigmasq,
       main = expression(sigma^2), xlab = "NNGP_Stan", ylab = "NNGP_fast")
abline(a = 0, b = 1)
qqplot(sim_w_post$beta[, 1], pos.beta[1, ],
       main = expression(beta[1]), xlab = "NNGP_Stan", ylab = "NNGP_fast")
abline(a = 0, b = 1)
qqplot(sim_w_post$beta[, 2], pos.beta[2, ],
       main = expression(beta[2]), xlab = "NNGP_Stan", ylab = "NNGP_fast")
abline(a = 0, b = 1)


par(mfrow = c(1, 1))
pdf(paste("./pic/qq_pseodu_RE_beta1.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")
qqplot(pos.beta[1, ], sim_w_post$beta[, 1]
       , xlab = "", ylab = "NNGP latent", main = "")
abline(a = 0, b = 1)
dev.off()

pdf(paste("./pic/qq_pseodu_RE_beta2.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")
qqplot( pos.beta[2, ], sim_w_post$beta[, 2]
        , xlab = "", ylab = "NNGP latent", main = "")
abline(a = 0, b = 1)
dev.off()

pdf(paste("./pic/qq_pseodu_RE_sigmasq.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")
qqplot(pos.sigmasq, sim_w_post$sigmasq
       , xlab = "", ylab = "NNGP latent", main = "")
abline(a = 0, b = 1)
dev.off()



par(mfrow = c(1, 1))
pdf(paste("./pic/qq_pseodu_RE_w1.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")
qqplot(pos.w[ 1, ], sim_w_post$w[, 1] - sim_w_post$beta[, 1]
       , xlab = "", ylab = "NNGP latent", main = "")
abline(a = 0, b = 1)
dev.off()

pdf(paste("./pic/qq_pseodu_RE_w2.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")
qqplot( pos.w[ 2, ], sim_w_post$w[, 2] - sim_w_post$beta[, 1]
       , xlab = "", ylab = "NNGP latent", main = "")
abline(a = 0, b = 1)
dev.off()

pdf(paste("./pic/qq_pseodu_RE_w3.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")
qqplot(pos.w[ 3, ], sim_w_post$w[, 3] - sim_w_post$beta[, 1]
       , xlab = "", ylab = "NNGP latent", main = "")
abline(a = 0, b = 1)
dev.off()


par(mfrow = c(1, 3))

qqplot(sim_w_post$sigma^2, pos.p.sigmasq,
       main = expression(sigma^2), xlab = "NNGP_Stan", ylab = "NNGP_fast_point")
abline(a = 0, b = 1)
qqplot(sim_w_post$beta[, 1], pos.p.beta[1, ],
       main = expression(beta[1]), xlab = "NNGP_Stan", ylab = "NNGP_fast_point")
abline(a = 0, b = 1)
qqplot(sim_w_post$beta[, 2], pos.p.beta[2, ],
       main = expression(beta[2]), xlab = "NNGP_Stan", ylab = "NNGP_fast_point")
abline(a = 0, b = 1)

par(mfrow = c(1, 1))
pdf(paste("./pic/qq_conj_RE_beta1.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")
qqplot(pos.p.beta[1, ], sim_w_post$beta[, 1]
       , xlab = "", ylab = "NNGP latent", main = "")
abline(a = 0, b = 1)
dev.off()


pdf(paste("./pic/qq_conj_RE_beta2.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")
qqplot(pos.p.beta[2, ], sim_w_post$beta[, 2]
       , xlab = "", ylab = "NNGP latent", main = ""
)
abline(a = 0, b = 1)
dev.off()

pdf(paste("./pic/qq_conj_RE_sigmasq.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")
qqplot(pos.p.sigmasq, sim_w_post$sigmasq
       , xlab = "", ylab = "NNGP latent", main = "" )
abline(a = 0, b = 1)
dev.off()



par(mfrow = c(1, 1))
pdf(paste("./pic/qq_conj_RE_w1.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")
qqplot(pos.p.w[ 1, ], sim_w_post$w[, 1] - sim_w_post$beta[, 1]
       , xlab = "", ylab = "NNGP latent", main = "")
abline(a = 0, b = 1)
dev.off()


pdf(paste("./pic/qq_conj_RE_w2.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")
qqplot(pos.p.w[ 2, ], sim_w_post$w[, 2] - sim_w_post$beta[, 1]
       , xlab = "", ylab = "NNGP latent", main = ""
       )
abline(a = 0, b = 1)
dev.off()

pdf(paste("./pic/qq_conj_RE_w3.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")
qqplot(pos.p.w[ 3, ], sim_w_post$w[, 3] - sim_w_post$beta[, 1]
       , xlab = "", ylab = "NNGP latent", main = "" )
abline(a = 0, b = 1)
dev.off()


#-------------------------summary of conj and pseudo-------------------------#
quantile(pos.p.beta[1, ], myqnt); mean(pos.p.beta[1, ])
quantile(pos.p.beta[2, ], myqnt); mean(pos.p.beta[2, ])
quantile(pos.p.sigmasq, myqnt); mean(pos.p.sigmasq)
quantile(pos.p.sigmasq, myqnt)*0.06; mean(pos.p.sigmasq*0.06)


quantile(pos.beta[1, ], myqnt); mean(pos.beta[1, ])
quantile(pos.beta[2, ], myqnt); mean(pos.beta[2, ])
quantile(pos.sigmasq, myqnt); mean(pos.sigmasq)
















































#title('Model1 (M = 15)', outer = TRUE, cex.main = 2, font.main = 4, line = -8)

#--------------------- MCMC chain comparation (phi = 30) ---------------------#

niter = 28000
### read output ###
library(coda)
returnSamples_fullGP <- matrix(
  (read.table("../RWM_version_model1/results/returnSamples_fullGP_unif.txt", 
              quote="\"", comment.char="")$V1), nrow = niter, ncol = 5, byrow = T)

returnSamples_nngp <- matrix(
  (read.table("../RWM_version_model1/results/returnSamples_nngp_unif.txt", 
              quote="\"", comment.char="")$V1), nrow = niter, ncol = 5, byrow = T)

niter2 = 300000
returnSamples_w_nngp <- matrix(
  (read.table(
    "./RWM_version_model2_flat/results_seed1234/sorted_coord_order/returnSamples_nngp_w.txt", 
    quote="\"", comment.char="")$V1), nrow = niter2, ncol = 5, byrow = T)


ap = 3; bp = 300
sigmasq_sample_fullGP <- returnSamples_fullGP[,4]^2
tausq_sample_fullGP <- returnSamples_fullGP[,3]^2
phi_sample_fullGP <- (bp + ap * exp(returnSamples_fullGP[, 5]))/
  (exp(returnSamples_fullGP[, 5]) + 1)
sigmasq_sample_nngp <- returnSamples_nngp[,4]^2
tausq_sample_nngp <- returnSamples_nngp[,3]^2
phi_sample_nngp <- (bp + ap * exp(returnSamples_nngp[, 5]))/
  (exp(returnSamples_nngp[, 5]) + 1)

phi_sample_w_nngp <- (bp + ap * exp(returnSamples_w_nngp[, 5]))/
  (exp(returnSamples_w_nngp[, 5]) + 1)


library("ggplot2")
library("gridExtra")
library("TeachingDemos")  #subplot
library("cowplot")
library("grid")

fullGP_data <- data.frame(ind = 1: 28000, beta1 = returnSamples_fullGP[, 1],
                        beta2 = returnSamples_fullGP[, 2], 
                        sigmasq = returnSamples_fullGP[,4]^2,
                        tausq = returnSamples_fullGP[,3]^2, 
                        phi = phi_sample_fullGP)

nngp_data <- data.frame(ind = 1: 28000, beta1 = returnSamples_nngp[, 1],
                        beta2 = returnSamples_nngp[, 2], 
                        sigmasq = returnSamples_nngp[,4]^2,
                        tausq = returnSamples_nngp[,3]^2, 
                        phi = phi_sample_nngp)

hmc_data <- data.frame(ind = 1: 2000, beta1 = sim_post_30[, 1, "beta[1]"],
                       beta2 = sim_post_30[, 1, "beta[2]"], 
                       sigmasq = sim_post_30[, 1, "sigma"]^2,
                       tausq = sim_post_30[, 1, "tau"]^2,
                       phi = sim_post_30[, 1, "phi"])


parameters = c("beta1", "beta2", "sigmasq", "tausq", "phi")
sparse_grid <- seq(0, 28000, by = 14)[-1]
sparse_grid2 <- seq(20000,22000, by = 14)
gender_colors <- c("mediumseagreen", "steelblue3","red")
names(gender_colors) <- c("fullGP_M-H", "nngp_M-H", "nngp_HMC")
for (i in 1: length(parameters)){
  temp_p <- parameters[i]
  p2 <- ggplot() + 
    #geom_line(data = fullGP_data[sparse_grid2, ], 
    #          aes(x = ind/20, y = eval(parse(text = temp_p))),
    #          color = "mediumseagreen", 
    #          linetype = 1, size = 0.8) +
    geom_line(data = hmc_data[sparse_grid2/14,  ], 
              aes(x = ind , y = eval(parse(text = temp_p))), color = "red", 
              linetype = 1, size = 0.8) +
    geom_line(data = nngp_data[sparse_grid2, ], 
              aes(x = ind/14, y = eval(parse(text = temp_p))),
              color = "steelblue3", 
              linetype = 1, size = 0.8)  +
    theme(axis.title.x=element_blank(), 
          axis.title.y=element_blank()) 
  p_s = ggplotGrob(p2)
  v <- range(nngp_data[, parameters[i]])
  v <- c(v[1] + (0.25) * (v[2] - v[1]),  v[2]-0.5)
  #v <- c(sum(v)/2, v[2])
  #v <- c(v[1], v[1] + (0.75) * (v[2] - v[1]))
  p1 <- ggplot(data = nngp_data[sparse_grid, ],
               aes(x = ind / 14, y = eval(parse(text = temp_p)))) +
    #geom_line(data = fullGP_data[sparse_grid,], 
    #          aes(x = ind/20, y = eval(parse(text = temp_p)), 
    #              color = "fullGP_M-H"), size = 0.2) +
    geom_line(data = hmc_data[sparse_grid/14, ], 
              aes(x = ind , y = eval(parse(text = temp_p)), 
                  color = "nngp_HMC"), size = 0.2) +
    geom_line(data = nngp_data[sparse_grid,],
              aes(x = ind/14, y = eval(parse(text = temp_p)), 
                  color = "nngp_M-H"), size = 0.2) +
    scale_color_manual(name="", values=gender_colors) +
    theme_bw() +
    theme(axis.title.x=element_blank(), 
          axis.title.y=element_blank(),
          plot.title = element_text(
            lineheight=.8, face="bold", size = 15)) +
    annotation_custom(grob = ggplotGrob(p2),
                      xmin = 200, xmax = 000, 
                      ymin = v[1], ymax = v[2]) + 
    ggtitle(temp_p) +xlab("ind")
  p1
  readline(prompt = "Pause. Press <Enter> to continue...")
}

















load("./stan_test_paper_flat/results_30_seed1234/sorted_coord_order/model2_nngp_10_6000")
samples_w_30 <- samples_w
sim_w_post_30 <- extract (samples_w_30, inc_warmup = T, permuted = F)


ind2 <- c(1: 300000)
nngp_data <- data.frame(ind = 1: 300000, beta1 = returnSamples_w_nngp[ind2, 1],
                        beta2 = returnSamples_w_nngp[ind2, 2], 
                        sigmasq = returnSamples_w_nngp[ind2, 4]^2,
                        tausq = returnSamples_w_nngp[ind2, 3]^2, 
                        phi = phi_sample_w_nngp[ind2])

hmc_data <- data.frame(ind = 1: 6000, beta1 = sim_w_post_30[, 1, "beta[1]"],
                       beta2 = sim_w_post_30[, 1, "beta[2]"], 
                       sigmasq = sim_w_post_30[, 1, "sigma"]^2,
                       tausq = sim_w_post_30[, 1, "tau"]^2,
                       phi = sim_w_post_30[, 1, "phi"])


jump = 43
parameters = c("beta1", "beta2", "sigmasq", "tausq", "phi")
sparse_grid <- seq(jump*35, jump*6000, by = jump)[-1]
sparse_grid2 <- seq(200000, 240000, by = jump)
gender_colors <- c("mediumseagreen", "steelblue3","red")
names(gender_colors) <- c("fullGP_M-H", "nngp_M-H", "nngp_HMC")
par(mfrow = c(1, 1))
for (i in 1: length(parameters)){
  temp_p <- parameters[i]
  p2 <- ggplot() + 
    #geom_line(data = fullGP_data[sparse_grid2, ], 
    #          aes(x = ind/20, y = eval(parse(text = temp_p))),
    #          color = "mediumseagreen", 
    #          linetype = 1, size = 0.8) +
    geom_line(data = hmc_data[sparse_grid2/jump , ], 
              aes(x = ind , y = eval(parse(text = temp_p))), color = "red", 
              linetype = 1, size = 0.8) +
    geom_line(data = nngp_data[sparse_grid2, ], 
              aes(x = ind/jump , y = eval(parse(text = temp_p))),
              color = "steelblue3", 
              linetype = 1, size = 0.8)  +
    theme(axis.title.x=element_blank(), 
          axis.title.y=element_blank()) 
  p_s = ggplotGrob(p2)
  #v <- range(nngp_data[sparse_grid, parameters[i]])
  v <- range(hmc_data[sparse_grid/jump , parameters[i]])
  v <- c(v[1] + (0.25) * (v[2] - v[1]), v[2])
  #v <- c(sum(v)/2, v[2])
  #v <- c(v[1], v[1] + (0.75) * (v[2] - v[1]))
  p1 <- ggplot(data = nngp_data[sparse_grid, ],
         aes(x = ind / jump , y = eval(parse(text = temp_p)))) +
    #ggplot(data = fullGP_data[sparse_grid, ],
    #           aes(x = ind / 50, y = eval(parse(text = temp_p)))) +
    #geom_line(data = fullGP_data[sparse_grid,], 
    #          aes(x = ind/20, y = eval(parse(text = temp_p)), 
    #              color = "fullGP_M-H"), size = 0.2) +
    geom_line(data = hmc_data[sparse_grid/jump, ], 
              aes(x = ind, y = eval(parse(text = temp_p)), 
                  color = "nngp_HMC"), size = 0.2) +
    geom_line(data = nngp_data[sparse_grid,],
              aes(x = ind / jump , y = eval(parse(text = temp_p)), 
                  color = "nngp_M-H"), size = 0.2) +
    scale_color_manual(name="", values=gender_colors) +
    theme_bw() +
    theme(axis.title.x=element_blank(), 
          axis.title.y=element_blank(),
          plot.title = element_text(
            lineheight=.8, face="bold", size = 15)) +
    annotation_custom(grob = ggplotGrob(p2),
                      xmin = 1000, xmax = 6000, 
                      ymin = v[1], ymax = v[2] ) + 
    ggtitle(temp_p) +xlab("ind")
  p1
  readline(prompt = "Pause. Press <Enter> to continue...")
}







jump = 43
parameters = c("beta1", "beta2", "sigmasq", "tausq", "phi")
sparse_grid1 <- seq(jump, jump*6000, by = jump)[-1]
sparse_grid <- seq(jump*34, jump*6000, by = jump)[-1]
sparse_grid2 <- seq(200000, 200000+200*jump, by = jump)
gender_colors <- c("mediumseagreen", "steelblue3","red")
names(gender_colors) <- c("fullGP_M-H", "nngp_M-H (every 43 samples)", "nngp_HMC  ")
temp_p <- parameters[1]
p2 <- ggplot() + 
  #geom_line(data = fullGP_data[sparse_grid2, ], 
  #          aes(x = ind/20, y = eval(parse(text = temp_p))),
  #          color = "mediumseagreen", 
  #          linetype = 1, size = 0.8) +
  geom_line(data = hmc_data[sparse_grid2/jump , ], 
            aes(x = ind , y = eval(parse(text = temp_p))), color = "red", 
            linetype = 1, size = 0.8) +
  geom_line(data = nngp_data[sparse_grid2, ], 
            aes(x = ind/jump , y = eval(parse(text = temp_p))),
            color = "steelblue3", 
            linetype = 1, size = 0.8)  +
  theme(axis.title.x=element_blank(), 
        axis.title.y=element_blank()) 
p_s = ggplotGrob(p2)
v <- range(nngp_data[sparse_grid1, parameters[1]])
v <- c(v[1] + (0.3) * (v[2] - v[1]), v[2])
#v <- c(sum(v)/2, v[2])
#v <- c(v[1], v[1] + (0.75) * (v[2] - v[1]))
p1 <- ggplot(data = nngp_data[sparse_grid1, ],
             aes(x = ind / jump , y = eval(parse(text = temp_p)))) +
  #ggplot(data = fullGP_data[sparse_grid, ],
  #           aes(x = ind / 50, y = eval(parse(text = temp_p)))) +
  #geom_line(data = fullGP_data[sparse_grid,], 
  #          aes(x = ind/20, y = eval(parse(text = temp_p)), 
  #              color = "fullGP_M-H"), size = 0.2) +
  geom_line(data = hmc_data[sparse_grid1/jump, ], 
            aes(x = ind, y = eval(parse(text = temp_p)), 
                color = "nngp_HMC  "), size = 0.2) +
  geom_line(data = nngp_data[sparse_grid1,],
            aes(x = ind / jump , y = eval(parse(text = temp_p)), 
                color = "nngp_M-H (every 43 samples)"), size = 0.2) +
  scale_color_manual(name="", values=gender_colors) +
  theme_bw() + theme(legend.position = "none") +
  theme(axis.title.x=element_blank(), 
        axis.title.y=element_blank(),
        plot.title = element_text(
          lineheight=.8, face="bold", size = 15)) +
  annotation_custom(grob = ggplotGrob(p2),
                    xmin = 1000, xmax = 6000, 
                    ymin = v[1], ymax = v[2] ) + 
  ggtitle(temp_p) +xlab("ind")




get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

temp_p2 <- parameters[4]
p22 <- ggplot() + 
  #geom_line(data = fullGP_data[sparse_grid2, ], 
  #          aes(x = ind/20, y = eval(parse(text = temp_p))),
  #          color = "mediumseagreen", 
  #          linetype = 1, size = 0.8) +
  geom_line(data = hmc_data[sparse_grid2/jump , ], 
            aes(x = ind , y = eval(parse(text = temp_p2))), color = "red", 
            linetype = 1, size = 0.8) +
  geom_line(data = nngp_data[sparse_grid2, ], 
            aes(x = ind/jump , y = eval(parse(text = temp_p2))),
            color = "steelblue3", 
            linetype = 1, size = 0.8)  +
  theme(axis.title.x=element_blank(), 
        axis.title.y=element_blank()) 
p_s2 = ggplotGrob(p22)
#v2 <- range(nngp_data[sparse_grid, parameters[4]])
v2 <- range(hmc_data[sparse_grid/jump , parameters[4]])
v2 <- c(v2[1] + (0.3) * (v2[2] - v2[1]), v2[2])
#v <- c(sum(v)/2, v[2])
#v <- c(v[1], v[1] + (0.75) * (v[2] - v[1]))
p12 <- ggplot(data = nngp_data[sparse_grid, ],
             aes(x = ind / jump , y = eval(parse(text = temp_p2)))) +
  #ggplot(data = fullGP_data[sparse_grid, ],
  #           aes(x = ind / 50, y = eval(parse(text = temp_p)))) +
  #geom_line(data = fullGP_data[sparse_grid,], 
  #          aes(x = ind/20, y = eval(parse(text = temp_p)), 
  #              color = "fullGP_M-H"), size = 0.2) +
  geom_line(data = hmc_data[sparse_grid/jump, ], 
            aes(x = ind, y = eval(parse(text = temp_p2)), 
                color = "nngp_HMC  "), size = 0.2) +
  geom_line(data = nngp_data[sparse_grid,],
            aes(x = ind / jump , y = eval(parse(text = temp_p2)), 
                color = "nngp_M-H (every 43 samples)"), size = 0.2) +
  scale_color_manual(name="", values=gender_colors) +
  theme_bw() +
  theme(axis.title.x=element_blank(), 
        axis.title.y=element_blank(),
        plot.title = element_text(
          lineheight=.8, face="bold", size = 15)) +
  annotation_custom(grob = ggplotGrob(p22),
                    xmin = 1000, xmax = 6000, 
                    ymin = v2[1], ymax = v2[2] ) + 
  ggtitle(temp_p2) +xlab("ind") + theme(legend.position = "top")

legend <- get_legend(p12)
p12 <- p12 + theme(legend.position="none")
grid.arrange(p1, p12, legend, ncol=2, nrow = 2,  
             layout_matrix = rbind(c(1,2), c(3,3)),
             widths = c(2.7, 2.7), heights = c(2.5, 0.2))

pdf("./stan_test_paper_flat/pic/HMCvsNNGP2.pdf", onefile = TRUE,
    width = 6, height =3.5 )
grid.arrange(p1, p12, legend, ncol=2, nrow = 2,  
             layout_matrix = rbind(c(1,2), c(3,3)),
             widths = c(4, 4), heights = c(2, 0.2))
dev.off()





