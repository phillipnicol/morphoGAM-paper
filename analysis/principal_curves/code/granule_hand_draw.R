
setwd(here::here("analysis/principal_curves/code"))

## Granule

xy <- readRDS("../../data/granule_slideseq.RDS")

## Prune outlier
xy.dist <- as.matrix(dist(xy))
nnk <- apply(xy.dist, 1, function(x) sort(x)[6])
outlier <- which(nnk > 2*median(nnk))
xy <- xy[-outlier,]

library(MorphoGAM)

ground_truth <- CurveFinderInteractive(xy, loop=FALSE)

save(ground_truth, file = "../data/granule_ground_truth.Rda")

#save curve and fitted coordinate
p <- ggarrange(ground_truth$curve.plot, 
          ground_truth$coordinate.plot, nrow=1, ncol=2, labels=c("a","b"))


ggsave(p, filename="../plots/granule_ground_truth.png", width=8, height=4)
