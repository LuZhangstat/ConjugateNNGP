# ---- draw pic ---- #
setwd("") # set to the path of ConjugateNNGP
setwd("asym_sim")
rm(list=ls())

library(rstan)

width <- 5
height <- 5
pointsize <- 16

#--- sim2 ---#
load("./results/sim2/model1_nngp_5")
load("./results/sim2/model2_nngp_5")
sam2M1m5 <- extract(samples)
sam2M2m5 <- extract(samples_w)

par(mfrow = c(1, 1))
pdf(paste("./pic/sim2m5qqphi.pdf", sep = ""),
    width=width, height=height, pointsize=pointsize, family="Courier")
qqplot(sam2M1m5$phi, sam2M2m5$phi, xlab = "response",
ylab = "latent", cex = 0.3)   # n_eff: 1291 vs 311
abline(a = 0, b = 1)
dev.off()

par(mfrow = c(1, 1))
pdf(paste("./pic/sim2m5qqdelta.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")
qqplot(sam2M1m5$tausq/sam2M1m5$sigmasq,   #n_eff: 945:1529 
       sam2M2m5b1$tausq/sam2M2m5$sigmasq, xlab = "response",
ylab = "latent", cex = 0.3)   #n_eff: 15:390
abline(a = 0, b = 1)
dev.off()


load("./results/sim2/model1_nngp_10")
load("./results/sim2/model2_nngp_10")
sam2M1m10 <- extract(samples)
sam2M2m10 <- extract(samples_w)

par(mfrow = c(1, 1))
pdf(paste("./pic/sim2m10qqphi.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")
qqplot(sam2M1m10$phi, sam2M2m10$phi, xlab = "response", 
ylab = "latent", cex = 0.3)   # n_eff: 1639 vs 402
abline(a = 0, b = 1)
dev.off()

pdf(paste("./pic/sim2m10qqdelta.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")
qqplot(sam2M1m10$tausq/sam2M1m10$sigmasq,   #n_eff: 970:1266
       sam2M2m10$tausq/sam2M2m10$sigmasq, xlab = "response", 
ylab = "latent", cex = 0.3)   #n_eff: 61:559
abline(a = 0, b = 1)
dev.off()



#--- sim3 ---#
load("./results/sim3/model1_nngp_5")
load("./results/sim3/model2_nngp_5")
sam3M1m5 <- extract(samples)
sam3M2m5 <- extract(samples_w)

par(mfrow = c(1, 1))
pdf(paste("./pic/sim3m5qqphi.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")
qqplot(sam3M1m5$phi, sam3M2m5$phi, xlab = "response", 
ylab = "latent", cex = 0.3)   # n_eff: 2178 vs 147
abline(a = 0, b = 1)
dev.off()

pdf(paste("./pic/sim3m5qqdelta.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")
qqplot(sam3M1m5$tausq/sam3M1m5$sigmasq,   #n_eff: 1943:2068 
       sam3M2m5$tausq/sam3M2m5$sigmasq, xlab = "response", 
ylab = "latent", cex = 0.3)   #n_eff: 47:237
abline(a = 0, b = 1)
dev.off()


load("./results/sim3/model1_nngp_10")
load("./results/sim3/model2_nngp_10")
sam3M1m10 <- extract(samples)
sam3M2m10 <- extract(samples_w)

par(mfrow = c(1, 1))
pdf(paste("./pic/sim3m10qqphi.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")
qqplot(sam3M1m10$phi, sam3M2m10$phi, xlab = "response", 
ylab = "latent", cex = 0.3)   # n_eff: 1737 vs 220
abline(a = 0, b = 1)
dev.off()

par(mfrow = c(1, 1))
pdf(paste("./pic/sim3m10qqdelta.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")
qqplot(sam3M1m10$tausq/sam3M1m10$sigmasq,   #n_eff: 1386:2105  
       sam3M2m10$tausq/sam3M2m10$sigmasq, xlab = "response", 
ylab = "latent", cex = 0.3)   #n_eff: 62:318
abline(a = 0, b = 1)
dev.off()

#--- sim4 ---#
load("./results/sim4/model1_nngp_5")
load("./results/sim4/model2_nngp_5")
sam4M1m5 <- extract(samples)
sam4M2m5 <- extract(samples_w)

par(mfrow = c(1, 1))
pdf(paste("./pic/sim4m5qqphi.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family = "Courier")
qqplot(sam4M1m5$phi, sam4M2m5$phi, xlab = "response",
ylab = "latent", cex = 0.3)
abline(a = 0, b = 1)
dev.off()

par(mfrow = c(1, 1))
pdf(paste("./pic/sim4m5qqdelta.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")
qqplot(sam4M1m5$tausq/sam4M1m5$sigmasq,
       sam4M2m5b1$tausq/sam4M2m5$sigmasq, xlab = "response",
ylab = "latent", cex = 0.3)
abline(a = 0, b = 1)
dev.off()


load("./results/sim4/model1_nngp_10")
load("./results/sim4/model2_nngp_10")
sam4M1m10 <- extract(samples)
sam4M2m10 <- extract(samples_w)

par(mfrow = c(1, 1))
pdf(paste("./pic/sim4m10qqphi.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")
qqplot(sam4M1m10$phi, sam4M2m10$phi, xlab = "response", 
ylab = "latent", cex = 0.3)   # n_eff: 1601 vs 203
abline(a = 0, b = 1)
dev.off()

par(mfrow = c(1, 1))
pdf(paste("./pic/sim4m10qqdelta.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")
qqplot(sam4M1m10$tausq/sam4M1m10$sigmasq,   #n_eff: 1544:1869
       sam4M2m10$tausq/sam4M2m10$sigmasq, xlab = "response", 
ylab = "latent", cex = 0.3)   #n_eff: 30:366
abline(a = 0, b = 1)
dev.off()


#--- sim5 ---#
load("./results/sim5/model1_nngp_5")
load("./results/sim5/model2_nngp_5")
load("./results/sim5/model2_nngp_5_b1")
sam5M1m5 <- extract(samples)
sam5M2m5 <- extract(samples_w)
sam5M2m5b1 <- extract(samples_w_b1)

par(mfrow = c(1, 1))
pdf(paste("./pic/sim5m5qqphi.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")
qqplot(sam5M1m5$phi, sam5M2m5$phi, xlab = "response", 
ylab = "latent", cex = 0.3)
abline(a = 0, b = 1)
dev.off()

par(mfrow = c(1, 1))
pdf(paste("./pic/sim5m5qqdelta.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")
qqplot(sam5M1m5$tausq/sam5M1m5$sigmasq,
       sam5M2m5$tausq/sam5M2m5$sigmasq, xlab = "response", 
ylab = "latent", cex = 0.3)
abline(a = 0, b = 1)
dev.off()

load("./results/sim5/model1_nngp_10")
load("./results/sim5/model2_nngp_10")
sam5M1m10 <- extract(samples)
sam5M2m10 <- extract(samples_w)

par(mfrow = c(1, 1))
pdf(paste("./pic/sim5m10qqphi.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")
qqplot(sam5M1m10$phi, sam5M2m10$phi, xlab = "response",
ylab = "latent", cex = 0.3)
abline(a = 0, b = 1)
dev.off()

par(mfrow = c(1, 1))
pdf(paste("./pic/sim5m10qqdelta.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")
qqplot(sam5M1m10$tausq/sam5M1m10$sigmasq,
       sam5M2m10b1$tausq/sam5M2m10$sigmasq, xlab = "response",
ylab = "latent", cex = 0.3)
abline(a = 0, b = 1)
dev.off()








