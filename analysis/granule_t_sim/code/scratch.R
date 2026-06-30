

#Load 
xy <- readRDS(here::here("analysis", "data", "granule_slideseq.RDS"))
xy <- as.matrix(xy)

xy.dist <- as.matrix(dist(xy))
nnk <- apply(xy.dist, 1, function(x) sort(x)[6])
outlier <- which(nnk > 2 * median(nnk))
if(length(outlier) > 0) {
  xy <- xy[-outlier, , drop = FALSE]
}

colnames(xy) <- c("xcoord", "ycoord")
fit <- CurveFinder(xy, knn = 9)

load("../../principal_curves/data/granule_ground_truth.Rda")
fit <- ground_truth

setwd(here::here("analysis", "granule_t_sim/code"))
Y_sim <- read.csv("../data/granule_simulation_peak/granule_simulated_dataset_peakt_1_sd.r_2_rep_05.csv",
row.names=1, check.names=FALSE)

Y_sim <- as.matrix(Y_sim)
storage.mode(Y_sim) <- "numeric"
rownames(Y_sim) <- as.character(seq_len(nrow(Y_sim)))



mgam <- MorphoGAM::MorphoGAM(Y_sim,
  curve.fit = fit,
  design = y ~ s(t) + s(r),
  return.fx = TRUE
)

library(tidyverse)
mgam$results |> arrange(desc(peak.t)) |> head()

library(tidyverse)
mgam$results |> arrange(desc(range.t)) |> head()

data.frame(x=fit$xyt$x, y=fit$xyt$y, color=Y_sim[1,]) |>
  ggplot(aes(x=x,y=y,color=color)) + geom_point() +
  scale_color_gradient(low="navyblue", high="firebrick1") +
  theme_bw() +
  ggtitle("MorphoGAM (knn=9, default)")
