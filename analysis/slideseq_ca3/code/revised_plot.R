setwd(here::here("analysis/slideseq_ca3/code"))

library(tidyverse)
library(mgcv)
library(ggpubr)
library(MorphoGAM)
library(scales)

library(STexampleData)

spe <- STexampleData::SlideSeqV2_mouseHPC()

ixs <- which(spe$celltype == "CA3")

xy <- spatialCoords(spe)[ixs,]
Y <- counts(spe)[,ixs]

xy.dist <- as.matrix(dist(xy))
knn <- 20
prune.outlier <- 3
nnk <- apply(xy.dist, 1, function(x) sort(x)[knn + 1])
outlier <- which(nnk > (prune.outlier) * median(nnk))

xy <- xy[-outlier,]
Y <- Y[, -outlier]

fit <- CurveFinder(xy)

load("../data/mgam_ca3_t_and_r.RData")
mgam <- mgam2

Y <- as.matrix(Y)

p.prev <- plotGAMestimates(
  Y,
  genes = c("Rgs14", "Cpne9"),
  mgam_object = mgam,
  curve_fit = fit
) +
  scale_y_sqrt()

p.peak <- plotGAMestimates(
  Y,
  genes = order(mgam$results$peak.t, decreasing = TRUE)[1:5],
  mgam_object = mgam,
  curve_fit = fit
) +
  scale_y_sqrt()

p.range <- plotGAMestimates(
  Y,
  genes = order(mgam$results$range.t, decreasing = TRUE)[1:5],
  mgam_object = mgam,
  curve_fit = fit
) +
  scale_y_sqrt()


### Pvalues and peak statistic
mgam$results["Rgs14",]
mgam$results["Cpne9",]
mgam$results["Fxyd6",]
mgam$results["Hpca",]

p.curve <- fit$curve.plot + ggtitle("CA3 cells")

sim.data.dir <- here::here("analysis", "ca3_newsim", "data")

sim.results.file <- file.path(sim.data.dir, "ca3_rmodel_results.csv")
if(!file.exists(sim.results.file)) {
  sim.results.file <- file.path(sim.data.dir, "ca3_rmodel_results_smoke.csv")
}

if(!file.exists(sim.results.file)) {
  stop("No CA3 simulation result file found in ", sim.data.dir)
}

sim.results <- read_csv(sim.results.file, show_col_types = FALSE)

spateo.results.file <- file.path(sim.data.dir, "ca3_spateo_model_results.csv")
if(file.exists(spateo.results.file)) {
  spateo.results <- read_csv(spateo.results.file, show_col_types = FALSE) |>
    select(pattern, peak.t.value, sd.r.value, replicate, spateo_rank)

  sim.results <- sim.results |>
    left_join(
      spateo.results,
      by = c("pattern", "peak.t.value", "sd.r.value", "replicate")
    )
}

gaston.results.file <- file.path(sim.data.dir, "ca3_gaston_model_results.csv")
if(file.exists(gaston.results.file)) {
  gaston.results <- read_csv(gaston.results.file, show_col_types = FALSE) |>
    select(pattern, peak.t.value, sd.r.value, replicate, gaston_rank)

  sim.results <- sim.results |>
    left_join(
      gaston.results,
      by = c("pattern", "peak.t.value", "sd.r.value", "replicate")
    )
}

sim.xy <- read_csv(file.path(sim.data.dir, "ca3_xy.csv"), show_col_types = FALSE)
sim.t <- read_csv(file.path(sim.data.dir, "ca3_t_coords.csv"), show_col_types = FALSE)[[1]]
sim.r <- read_csv(file.path(sim.data.dir, "ca3_r_coords.csv"), show_col_types = FALSE)[[1]]

sim.rank.cols <- c(
  "mgam_rank",
  "spark_rank",
  "nnsvg_rank",
  "heartsvg_rank",
  "spateo_rank",
  "gaston_rank"
)

missing.rank.cols <- setdiff(sim.rank.cols, colnames(sim.results))
if(length(missing.rank.cols) > 0) {
  stop("Missing rank columns in ", sim.results.file, ": ", paste(missing.rank.cols, collapse = ", "))
}

sim.method.labels <- c(
  mgam_rank = "MorphoGAM",
  spark_rank = "SPARK-X",
  nnsvg_rank = "nnSVG",
  heartsvg_rank = "HEARTSVG",
  spateo_rank = "Spateo",
  gaston_rank = "GASTON"
)

sim.method.palette <- c(
  "MorphoGAM" = "#0072B2",
  "SPARK-X" = "#009E73",
  "nnSVG" = "#D55E00",
  "HEARTSVG" = "#CC79A7",
  "Spateo" = "#E69F00",
  "GASTON" = "#56B4E9"
)

rescale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if(diff(rng) == 0) {
    return(rep(0, length(x)))
  }
  (x - rng[1]) / diff(rng)
}

localized.signal <- rescale01(exp(-0.5 * ((sim.t - 0.5) / 0.12)^2))

sim.schematic <- tibble(
  xcoord = sim.xy$xcoord,
  ycoord = sim.xy$ycoord,
  signal = localized.signal,
  confounder = rescale01(sim.r)
)

spatial.theme <- theme_bw(base_size = 11) +
  theme(
    axis.title = element_text(),
    axis.text = element_text(color = "black"),
    legend.position = "right",
    legend.title = element_text(),
    plot.title = element_text(hjust = 0),
    plot.margin = margin(4, 4, 4, 4)
  )

make_spatial_panel <- function(df, value_col, title, legend_title) {
  ggplot(df, aes(x = xcoord, y = ycoord, color = .data[[value_col]])) +
    geom_point(size = 0.6, alpha = 0.95) +
    coord_equal() +
    scale_color_viridis_c(option = "magma", name = legend_title) +
    labs(title = title, x = "x", y = "y") +
    guides(color = guide_colorbar(barheight = unit(26, "pt"))) +
    spatial.theme
}

p.signal <- make_spatial_panel(
  sim.schematic,
  value_col = "signal",
  title = "Localized signal",
  legend_title = "Signal"
)

p.confounder <- make_spatial_panel(
  sim.schematic,
  value_col = "confounder",
  title = "Confounding genes 2 - G: linear radial signal",
  legend_title = "Signal"
)

sim.rank.df <- sim.results |>
  filter(.data$pattern == "peak") |>
  select(pattern, peak.t.value, sd.r.value, replicate, all_of(sim.rank.cols)) |>
  pivot_longer(
    cols = all_of(sim.rank.cols),
    names_to = "method",
    values_to = "rank"
  ) |>
  mutate(
    method = recode(.data$method, !!!sim.method.labels),
    method = factor(.data$method, levels = unname(sim.method.labels)),
    sd.r.value = factor(
      .data$sd.r.value,
      levels = c(0.5, 1, 2),
      labels = c("Low r", "Medium r", "High r")
    )
  ) |>
  filter(!is.na(.data$rank), .data$rank > 0)

p.rank <- ggplot(sim.rank.df, aes(x = method, y = rank, fill = method)) +
  geom_boxplot(width = 0.62, outlier.shape = NA, linewidth = 0.35, alpha = 0.88) +
  geom_point(
    aes(color = method),
    position = position_jitter(width = 0.12, height = 0, seed = 1),
    size = 1.45,
    alpha = 0.85
  ) +
  facet_wrap(~ sd.r.value, nrow = 1) +
  scale_fill_manual(values = sim.method.palette) +
  scale_color_manual(values = sim.method.palette) +
  scale_y_sqrt(
    breaks = pretty_breaks(n = 5),
    expand = expansion(mult = c(0.05, 0.12))
  ) +
  labs(x = NULL, y = "Rank of gene 1") +
  theme_bw(base_size = 11) +
  theme(
    axis.title = element_text(),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 35, hjust = 1),
    strip.background = element_rect(fill = "grey92", color = "grey70"),
    strip.text = element_text(),
    legend.position = "none",
    plot.margin = margin(4, 4, 4, 4),
    panel.spacing.x = unit(0.9, "lines")
  )

p.localized.sim <- ggarrange(
  ggarrange(p.signal, p.confounder, nrow = 1, ncol = 2),
  p.rank,
  nrow = 2,
  heights = c(1, 1.15),
  labels = c("d", "e")
)

p.full.revised <- ggarrange(
  ggarrange(
    p.curve,
    p.prev,
    nrow = 1,
    ncol = 2,
    labels = c("a", "b", "")
  ),
  ggarrange(
    p.peak,
    p.range,
    nrow = 2
  ),
  p.localized.sim,
  nrow = 3,
  heights = c(1.5, 2, 2.7),
  labels = c("", "c", "")
)

ggsave(
  p.full.revised,
  filename = "../plots/ca3_full_plot_revised.png",
  width = 9.11 * 1.1,
  height = 13.8,
  units = "in"
)

ggsave(
  p.full.revised,
  filename = "../plots/ca3_full_plot_revised.pdf",
  width = 9.11 * 1.1,
  height = 13.8,
  units = "in"
)
