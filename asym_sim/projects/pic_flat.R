# ---- draw pic ---- #
setwd("h:/research/proj1/asym_sim")
rm(list=ls())

library(rstan)
# library(rstanarm)

width <- 5
height <- 5
pointsize <- 16



#--- sim1 ---#
load("./results/sim1/model1_nngp_5")
load("./results/sim1/model2_nngp_5")
load("./results/sim1/model2_nngp_5_b1")
sam1M1m5 <- extract(samples)
sam1M2m5 <- extract(samples_w)
sam1M2m5b1 <- extract(samples_w_b1)


par(mfrow = c(1, 1))
pdf(paste("./pic/sim1m5qqphi.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")
qqplot(sam1M1m5$phi, sam1M2m5$phi, xlab = "response", 
ylab = "latent", cex = 0.3)   # n_eff: 723 vs 102
abline(a = 0, b = 1)
dev.off()

par(mfrow = c(1, 1))
pdf(paste("./pic/sim1m5qqdelta.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")
qqplot(sam1M1m5$tausq/sam1M1m5$sigmasq,   #n_eff: 398:492 
       sam1M2m5$tausq/sam1M2m5$sigmasq, xlab = "response", 
ylab = "latent", cex = 0.3)   #n_eff: 43:40
abline(a = 0, b = 1)
dev.off()


#qqplot(sam1M1m5$phi, sam1M2m5b1$phi)   # n_eff: 723 vs 836
#abline(a = 0, b = 1)

#qqplot(log(sam1M1m5$tausq/sam1M1m5$sigmasq),   #n_eff: 398:492 
#       log(sam1M2m5b1$tausq/sam1M2m5b1$sigmasq))   #n_eff: 24:23
#abline(a = 0, b = 1)




load("./results/sim1/model1_nngp_10")
load("./results/sim1/model2_nngp_10")
load("./results/sim1/model2_nngp_10_b1")
sam1M1m10 <- extract(samples)
sam1M2m10 <- extract(samples_w)
sam1M2m10b1 <- extract(samples_w_b1)

par(mfrow = c(1, 2))

#qqplot(sam1M1m10$phi, sam1M2m10$phi)   # n_eff: 1082 vs 71
#abline(a = 0, b = 1)

#qqplot(log(sam1M1m10$tausq/sam1M1m10$sigmasq),   #n_eff: 413:486 
#       log(sam1M2m10$tausq/sam1M2m10$sigmasq))   #n_eff: 28:13
#abline(a = 0, b = 1)


par(mfrow = c(1, 1))
pdf(paste("./pic/sim1m10qqphi.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")
qqplot(sam1M1m10$phi, sam1M2m10b1$phi, xlab = "response", 
ylab = "latent", cex = 0.3)   # n_eff: 1082 vs 472
abline(a = 0, b = 1)
dev.off()

par(mfrow = c(1, 1))
pdf(paste("./pic/sim1m10qqdelta.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")
qqplot(sam1M1m10$tausq/sam1M1m10$sigmasq,   #n_eff: 413:486 
       sam1M2m10b1$tausq/sam1M2m10b1$sigmasq, xlab = "response", 
ylab = "latent", cex = 0.3)   #n_eff: 38:38
abline(a = 0, b = 1)
dev.off()


#--- sim2 ---#
load("./results/sim2/model1_nngp_5")
load("./results/sim2/model2_nngp_5")
load("./results/sim2/model2_nngp_5_b1")
sam2M1m5 <- extract(samples)
sam2M2m5 <- extract(samples_w)
sam2M2m5b1 <- extract(samples_w_b1)

par(mfrow = c(1, 2))

#qqplot(sam2M1m5$phi, sam2M2m5$phi)   # n_eff: 1291 vs 85
#abline(a = 0, b = 1)

#qqplot(sam2M1m5$tausq/sam2M1m5$sigmasq,   #n_eff: 945:1529 
#       sam2M2m5$tausq/sam2M2m5$sigmasq)   #n_eff: 36:29
#abline(a = 0, b = 1)

par(mfrow = c(1, 1))
pdf(paste("./pic/sim2m5qqphi.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")
qqplot(sam2M1m5$phi, sam2M2m5b1$phi, xlab = "response", 
ylab = "latent", cex = 0.3)   # n_eff: 1291 vs 311
abline(a = 0, b = 1)
dev.off()

par(mfrow = c(1, 1))
pdf(paste("./pic/sim2m5qqdelta.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")
qqplot(sam2M1m5$tausq/sam2M1m5$sigmasq,   #n_eff: 945:1529 
       sam2M2m5b1$tausq/sam2M2m5b1$sigmasq, xlab = "response", 
ylab = "latent", cex = 0.3)   #n_eff: 15:390
abline(a = 0, b = 1)
dev.off()


load("./results/sim2/model1_nngp_10")
load("./results/sim2/model2_nngp_10")
load("./results/sim2/model2_nngp_10_b1")
sam2M1m10 <- extract(samples)
sam2M2m10 <- extract(samples_w)
sam2M2m10b1 <- extract(samples_w_b1)

par(mfrow = c(1, 2))
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

#qqplot(sam2M1m10$phi, sam2M2m10b1$phi)   # n_eff: 1639 vs 283
#abline(a = 0, b = 1)

#qqplot(log(sam2M1m10$tausq/sam2M1m10$sigmasq),   #n_eff: 970:1266
#       log(sam2M2m10b1$tausq/sam2M2m10b1$sigmasq))   #n_eff: 39:482
#abline(a = 0, b = 1)


#--- sim3 ---#
load("./results/sim3/model1_nngp_5")
load("./results/sim3/model2_nngp_5")
load("./results/sim3/model2_nngp_5_b1")
sam3M1m5 <- extract(samples)
sam3M2m5 <- extract(samples_w)
sam3M2m5b1 <- extract(samples_w_b1)

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

qqplot(sam3M1m5$phi, sam3M2m5b1$phi, cex = 0.3)   # n_eff: 2178 vs 431
abline(a = 0, b = 1)

#qqplot(sam3M1m5$tausq/sam3M1m5$sigmasq,   #n_eff: 1943:2068 
#       sam3M2m5b1$tausq/sam3M2m5b1$sigmasq, cex = 0.3)   #n_eff: 87:3000
#abline(a = 0, b = 1)


load("./results/sim3/model1_nngp_10")
load("./results/sim3/model2_nngp_10")
load("./results/sim3/model2_nngp_10_b1")
sam3M1m10 <- extract(samples)
sam3M2m10 <- extract(samples_w)
sam3M2m10b1 <- extract(samples_w_b1)

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

#qqplot(sam3M1m10$phi, sam3M2m10b1$phi)   # n_eff: 1737 vs 239
#abline(a = 0, b = 1)

#qqplot(sam3M1m10$tausq/sam3M1m10$sigmasq,   #n_eff: 1386:2105  
#       sam3M2m10b1$tausq/sam3M2m10b1$sigmasq)   #n_eff: 67:320
#abline(a = 0, b = 1)


#--- sim4 ---#
load("./results/sim4/model1_nngp_5")
load("./results/sim4/model2_nngp_5")
load("./results/sim4/model2_nngp_5_b1")
sam4M1m5 <- extract(samples)
sam4M2m5 <- extract(samples_w)
sam4M2m5b1 <- extract(samples_w_b1)

par(mfrow = c(1, 2))

#qqplot(sam4M1m5$phi, sam4M2m5$phi)   # n_eff: 2144 vs 169
#abline(a = 0, b = 1)

#qqplot(log(sam4M1m5$tausq/sam4M1m5$sigmasq),   #n_eff: 1833:1829 
#       log(sam4M2m5$tausq/sam4M2m5$sigmasq))   #n_eff: 23:182  
#abline(a = 0, b = 1)

par(mfrow = c(1, 1))
pdf(paste("./pic/sim4m5qqphi.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family = "Courier")
qqplot(sam4M1m5$phi, sam4M2m5b1$phi, xlab = "response", 
ylab = "latent", cex = 0.3)   # n_eff: 2144 vs 119
abline(a = 0, b = 1)
dev.off()

par(mfrow = c(1, 1))
pdf(paste("./pic/sim4m5qqdelta.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")
qqplot(sam4M1m5$tausq/sam4M1m5$sigmasq,   #n_eff: 1833:1829 
       sam4M2m5b1$tausq/sam4M2m5b1$sigmasq, xlab = "response", 
ylab = "latent", cex = 0.3)   #n_eff: 27:485  
abline(a = 0, b = 1)
dev.off()


load("./results/sim4/model1_nngp_10")
load("./results/sim4/model2_nngp_10")
load("./results/sim4/model2_nngp_10_b1")
sam4M1m10 <- extract(samples)
sam4M2m10 <- extract(samples_w)
sam4M2m10b1 <- extract(samples_w_b1)

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


#qqplot(sam4M1m10$phi, sam4M2m10b1$phi)   # n_eff: 1601 vs 38
#abline(a = 0, b = 1)

#qqplot(log(sam4M1m10$tausq/sam4M1m10$sigmasq),   #n_eff: 1544:1869
#       log(sam4M2m10b1$tausq/sam4M2m10b1$sigmasq))   #n_eff: 9:23
#abline(a = 0, b = 1)


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
ylab = "latent", cex = 0.3)   # n_eff: 2064 vs 129
abline(a = 0, b = 1)
dev.off()

par(mfrow = c(1, 1))
pdf(paste("./pic/sim5m5qqdelta.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")
qqplot(sam5M1m5$tausq/sam5M1m5$sigmasq,   #n_eff: 2357:2019 
       sam5M2m5$tausq/sam5M2m5$sigmasq, xlab = "response", 
ylab = "latent", cex = 0.3)   #n_eff: 38:794
abline(a = 0, b = 1)
dev.off()

#qqplot(sam5M1m5$phi, sam5M2m5b1$phi)   # n_eff: 2064 vs 196
#abline(a = 0, b = 1)

#qqplot(sam5M1m5$tausq/sam5M1m5$sigmasq,   #n_eff: 2357:2019 
#       sam5M2m5b1$tausq/sam5M2m5b1$sigmasq)   #n_eff: 59:874
#abline(a = 0, b = 1)



load("./results/sim5/model1_nngp_10")
load("./results/sim5/model2_nngp_10")
load("./results/sim5/model2_nngp_10_b1")
sam5M1m10 <- extract(samples)
sam5M2m10 <- extract(samples_w)
sam5M2m10b1 <- extract(samples_w_b1)

par(mfrow = c(1, 2))

#qqplot(sam5M1m10$phi, sam5M2m10$phi)   # n_eff: 2121 vs 135
#abline(a = 0, b = 1)

#qqplot(sam5M1m10$tausq/sam5M1m10$sigmasq,   #n_eff: 2089:2187
#       sam5M2m10$tausq/sam5M2m10$sigmasq)   #n_eff: 32:693
#abline(a = 0, b = 1)

par(mfrow = c(1, 1))
pdf(paste("./pic/sim5m10qqphi.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")
qqplot(sam5M1m10$phi, sam5M2m10b1$phi, xlab = "response", 
ylab = "latent", cex = 0.3)   # n_eff: 2121 vs 209
abline(a = 0, b = 1)
dev.off()

par(mfrow = c(1, 1))
pdf(paste("./pic/sim5m10qqdelta.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")
qqplot(sam5M1m10$tausq/sam5M1m10$sigmasq,   #n_eff: 2089:2187
       sam5M2m10b1$tausq/sam5M2m10b1$sigmasq, xlab = "response", 
ylab = "latent", cex = 0.3)   #n_eff: 51:3000
abline(a = 0, b = 1)
dev.off()


#--- sim6 ---#
load("./results/sim6/model1_nngp_5")
load("./results/sim6/model2_nngp_5")
load("./results/sim6/model2_nngp_5_b1")
sam6M1m5 <- extract(samples)
sam6M2m5 <- extract(samples_w)
sam6M2m5b1 <- extract(samples_w_b1)

par(mfrow = c(1, 1))
pdf(paste("./pic/sim6m5qqphi.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")
qqplot(sam6M1m5$phi, sam6M2m5$phi, xlab = "response", 
ylab = "latent", cex = 0.3)   # n_eff: 1865 vs 223
abline(a = 0, b = 1)
dev.off()

par(mfrow = c(1, 1))
pdf(paste("./pic/sim6m5qqdelta.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")
qqplot(sam6M1m5$tausq/sam6M1m5$sigmasq,   #n_eff: 2643:1999
       sam6M2m5$tausq/sam6M2m5$sigmasq, xlab = "response", 
ylab = "latent", cex = 0.3)   #n_eff: 105:3000
abline(a = 0, b = 1)
dev.off()


par(mfrow = c(1, 2))
qqplot(sam6M1m5$phi, sam6M2m5b1$phi)   # n_eff: 1865 vs 373
abline(a = 0, b = 1)

qqplot(sam6M1m5$tausq/sam6M1m5$sigmasq,   #n_eff: 2643:1999
       sam6M2m5b1$tausq/sam6M2m5b1$sigmasq)   #n_eff: 93:3000
abline(a = 0, b = 1)


load("./results/sim6/model1_nngp_10")
load("./results/sim6/model2_nngp_10")
load("./results/sim6/model2_nngp_10_b1")
sam6M1m10 <- extract(samples)
sam6M2m10 <- extract(samples_w)
sam6M2m10b1 <- extract(samples_w_b1)

par(mfrow = c(1, 2))
par(mfrow = c(1, 1))
pdf(paste("./pic/sim6m10qqphi.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")
qqplot(sam6M1m10$phi, sam6M2m10$phi, xlab = "response", 
ylab = "latent", cex = 0.3)   # n_eff: 1965 vs 471
abline(a = 0, b = 1)
dev.off()

par(mfrow = c(1, 1))
pdf(paste("./pic/sim6m10qqdelta.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")
qqplot(sam6M1m10$tausq/sam6M1m10$sigmasq,   #n_eff: 2541:2106
       sam6M2m10$tausq/sam6M2m10$sigmasq, xlab = "response", 
ylab = "latent", cex = 0.3)   #n_eff: 155:2234
abline(a = 0, b = 1)
dev.off()


#par(mfrow = c(1, 2))

#qqplot(sam6M1m10$phi, sam6M2m10b1$phi)   # n_eff: 1965 vs 517
#abline(a = 0, b = 1)

#qqplot(log(sam6M1m10$tausq/sam6M1m10$sigmasq),   #n_eff: 2541:2106 
#       log(sam6M2m10b1$tausq/sam6M2m10b1$sigmasq))   #n_eff: 145:3000
#abline(a = 0, b = 1)

print(samples_w, pars = c("sigmasq", "tausq", "phi"))

print(samples_w_b1, pars = c("sigmasq", "tausq", "phi"))
















# size for the mapping of w               
width <- 5
height <- 5
pointsize <- 16

par(mfrow = c(2, 2))

qqplot(sim_w_post$sigmasq, pos.sigmasq,
       main = expression(sigma^2), xlab = "NNGP_Stan", ylab = "NNGP_fast")
abline(a = 0, b = 1)
qqplot(sim_w_post$beta[, 1], pos.beta[1, ],
       main = expression(beta[1]), xlab = "NNGP_Stan", ylab = "NNGP_fast")
abline(a = 0, b = 1)
qqplot(sim_w_post$beta[, 2], pos.beta[2, ],
       main = expression(beta[2]), xlab = "NNGP_Stan", ylab = "NNGP_fast")
abline(a = 0, b = 1)

par(mfrow = c(1, 1))
pdf(paste("./pic/qq_pseodu_RE_w1.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")
qqplot(pos.w[ 1, ], sim_w_post$w[, 1] - sim_w_post$beta[, 1]
       , xlab = "", ylab = "NNGP random effect", main = "")
abline(a = 0, b = 1)
dev.off()

pdf(paste("./pic/qq_pseodu_RE_w2.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")
qqplot( pos.w[ 2, ], sim_w_post$w[, 2] - sim_w_post$beta[, 1]
       , xlab = "", ylab = "NNGP random effect", main = "")
abline(a = 0, b = 1)
dev.off()

pdf(paste("./pic/qq_pseodu_RE_w3.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")
qqplot(sim_w_post$w[, 3] - sim_w_post$beta[, 1], pos.w[ 3, ]
       , xlab = "", ylab = "NNGP random effect", main = "")
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
pdf(paste("stan_test_paper_flat/pic/qq_conj_RE_w1.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")
qqplot(pos.p.w[ 1, ], sim_w_post$w[, 1] - sim_w_post$beta[, 1]
       , xlab = "conjugate", ylab = "NNGP random effect", main = "")
abline(a = 0, b = 1)
dev.off()


pdf(paste("stan_test_paper_flat/pic/qq_conj_RE_w2.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")
qqplot(pos.p.w[ 2, ], sim_w_post$w[, 2] - sim_w_post$beta[, 1]
       , xlab = "conjugate", ylab = "NNGP random effect", main = ""
       )
abline(a = 0, b = 1)
dev.off()

pdf(paste("stan_test_paper_flat/pic/qq_conj_RE_w3.pdf", sep=""), 
    width=width, height=height, pointsize=pointsize, family="Courier")
qqplot(pos.p.w[ 3, ], sim_w_post$w[, 3] - sim_w_post$beta[, 1]
       , xlab = "conjugate", ylab = "NNGP random effect", main = "" )
abline(a = 0, b = 1)
dev.off()


























