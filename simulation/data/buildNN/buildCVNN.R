setwd("./") # set to the path of ConjugateNNGP
setwd("./simulation")
rm(list = ls())
library(MBA)
library(fields)
source("./src/NNMatrix.R")
load("./data/simdata/simdata.RData")

k.fold = 5; M = 10
ord <- order(coords[, 1]) # sort by longitude
coords <- coords[ord, ]


## Build NN matrixs for Cross Validation
set.seed(123)
require(caret)
flds <- createFolds(1:N, k = k.fold, list = TRUE, returnTrain = FALSE)

NN.matrix <- list()
nn.mod.ho <- list()

for (i in 1: k.fold) {
  cat("\n", "Building NN matrics for ", i, "th folder")
  t <- proc.time()
  ho.id <- flds[[i]]
  coords.ho <- coords[ho.id, ]
  coords.mod <- coords[-ho.id, ]

# Build NN index for training data in ith folder
  NN.matrix[[i]] <- NNMatrix(coords = coords.mod, n.neighbors = M, 
                             n.omp.threads = 2, search.type = "brute")

# Build nearest neighbor index for holdout locations
  nn.mod.ho[[i]] <- nn2(coords.mod, coords.ho, k = M)
  nn.mod.ho[[i]]$NN_distM <- 
    sapply(1:nrow(coords.ho), ho_dist, nn.mod.ho[[i]]$nn.idx, coords.mod)
  cat("\t", " takes time: ", proc.time() - t)
}

save(file = "./data/buildNN/CVNNmatric.RData", 
     list = c("M", "k.fold", "flds", "NN.matrix", "nn.mod.ho"))
