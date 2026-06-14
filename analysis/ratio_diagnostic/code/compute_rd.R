
#Compute ratio diagnostic for 3 real datasets and a square

setwd(here::here("analysis/ratio_diagnostic/code"))

# Granule 

library(MorphoGAM)

xy <- readRDS("../../data/granule_slideseq.RDS")

## Prune outlier
xy.dist <- as.matrix(dist(xy))
nnk <- apply(xy.dist, 1, function(x) sort(x)[6])
outlier <- which(nnk > 2*median(nnk))
xy <- xy[-outlier,]


fit.granule <- CurveFinder(xy,knn=5)
cat("Ratio diagonistic is ", fit.granule$t_v_r_span, "\n")

# Colon 


meta <- read.csv("../../data/GL2_distal_colon_cell_type_and_locations_2023.08.11.csv")

library(tidyverse)
meta <- as.data.frame(meta)
meta.sub <- meta |> filter(slice_full_name == "20220518_WT_dcol_slice_3") |>
  filter(spatial_neighborhood_v1 == "Mucosa") |>
  filter(leiden_combined_v2 == "Enterocyte")
locus <- as.matrix(meta.sub[,c("x","y")])
fit.colon <- CurveFinder(locus,knn=10,loop=TRUE)
cat("Ratio diagonistic is ", fit.colon$t_v_r_span, "\n")

#CA3 cells 
spe <- STexampleData::SlideSeqV2_mouseHPC()
spe <- STexampleData::SlideSeqV2_mouseHPC()

ixs <- which(spe$celltype == "CA3") #subset to CA3

xy <- spatialCoords(spe)[ixs,]

xy.dist <- as.matrix(dist(xy))
knn <- 20
prune.outlier <- 3
nnk <- apply(xy.dist, 1, function(x) sort(x)[knn+1])
outlier <- which(nnk > (prune.outlier)*median(nnk))

xy <- xy[-outlier,]
df <- xy

fit.ca3 <- CurveFinder(xy) #Fit the curve
cat("Ratio diagonistic is ", fit.ca3$t_v_r_span, "\n")

set.seed(1)
n <- 1000
#Generate points on a square 
x <- runif(n, 0, 1)
y <- runif(n, 0, 1)

xy <- cbind(x, y)

fit.square <- CurveFinder(xy)
cat("Ratio diagonistic is ", fit.square$t_v_r_span, "\n")

