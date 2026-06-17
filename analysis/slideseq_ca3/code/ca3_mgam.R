setwd(here::here("analysis/slideseq_ca3/code"))

library(dplyr)
library(STexampleData)
library(igraph)
library(tidyverse)
library(MorphoGAM)

spe <- STexampleData::SlideSeqV2_mouseHPC()

ixs <- which(spe$celltype == "CA3") #subset to CA3

xy <- spatialCoords(spe)[ixs,]
Y <- counts(spe)[,ixs]

xy.dist <- as.matrix(dist(xy))
knn <- 20
prune.outlier <- 3
nnk <- apply(xy.dist, 1, function(x) sort(x)[knn+1])
outlier <- which(nnk > (prune.outlier)*median(nnk))

xy <- xy[-outlier,]; Y <- Y[,-outlier]

fit <- CurveFinder(xy)

Y <- as.matrix(Y)
mgam <- MorphoGAM(Y, curve.fit=fit,
                  design=y~s(t,bs="cr"))

save(mgam, file="../data/mgam_ca3.RData")



mgam2 <- MorphoGAM(Y, curve.fit=fit,
                  design=y~s(t,bs="cr") + s(r, bs="cr"))

MorphoGAM::plotGAMestimates(Y, genes=c("Fxyd6", "Cpne9"), mgam_object = mgam2,
                             curve_fit = fit) + scale_y_sqrt()

save(mgam2, file="../data/mgam_ca3_t_and_r.RData")

#Top 1000 varaible genes
my.var <- apply(Y, 1, var)
Y.sub <- Y[order(my.var, decreasing = TRUE)[1:1000],]
dim(Y.sub)


mgam.small <- MorphoGAM(Y.sub, curve.fit=fit,
                  design=y~s(t, bs="cc") + s(r, bs="tp"))
