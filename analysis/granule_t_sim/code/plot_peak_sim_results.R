setwd(here::here("analysis/granule_t_sim/code"))

project.root <- normalizePath(file.path("..", "..", ".."), mustWork = TRUE)
r.version.short <- paste0(
  "R-",
  R.version$major,
  ".",
  strsplit(R.version$minor, ".", fixed = TRUE)[[1]][1]
)
renv.library <- file.path(project.root, "renv", "library", r.version.short, R.version$platform)
if(dir.exists(renv.library)) {
  .libPaths(c(renv.library, .libPaths()))
}

library(dplyr)
library(ggpubr)
library(ggplot2)
library(patchwork)
library(readr)
library(scales)
library(splines)
library(tidyr)

data.dir <- file.path("..", "data")
plot.dir <- file.path("..", "plots")
dir.create(plot.dir, showWarnings = FALSE, recursive = TRUE)

pick_existing_file <- function(candidates, label) {
  existing <- candidates[file.exists(candidates)]
  if(length(existing) == 0) {
    message("No ", label, " file found. Checked: ", paste(candidates, collapse = ", "))
    return(NULL)
  }

  existing[[1]]
}

read_rank_file <- function(file, rank_cols) {
  if(is.null(file)) {
    return(NULL)
  }

  message("Reading ", file)
  cols <- names(read_csv(file, n_max = 0, show_col_types = FALSE))
  keep <- intersect(c("pattern", "peak.t.value", "sd.r.value", "replicate", rank_cols), cols)

  read_csv(file, show_col_types = FALSE) %>%
    select(all_of(keep))
}

join_results <- function(...) {
  pieces <- list(...) %>%
    Filter(Negate(is.null), .)

  if(length(pieces) == 0) {
    stop("No simulation result files were found.")
  }

  Reduce(
    function(x, y) {
      full_join(
        x,
        y,
        by = c("pattern", "peak.t.value", "sd.r.value", "replicate")
      )
    },
    pieces
  )
}

r.results.file <- pick_existing_file(
  file.path(
    data.dir,
    c(
      "peak_only_3.csv",
      "granule_rmodel_results_peak.csv",
      "granule_rmodel_results.csv",
      "granule_rmodel_results_peak_smoke.csv",
      "granule_rmodel_results_smoke.csv"
    )
  ),
  "R model results"
)

spateo.results.file <- pick_existing_file(
  file.path(
    data.dir,
    c(
      "granule_spateo_model_results_peak.csv",
      "granule_spateo_model_results.csv"
    )
  ),
  "Spateo results"
)

gaston.results.file <- pick_existing_file(
  file.path(
    data.dir,
    c(
      "granule_gaston_model_results_peak.csv",
      "granule_gaston_model_results.csv"
    )
  ),
  "GASTON results"
)

r.results <- read_rank_file(
  r.results.file,
  c("mgam_rank", "spark_rank", "nnsvg_rank", "heartsvg_rank")
)
spateo.results <- read_rank_file(spateo.results.file, "spateo_rank")
gaston.results <- read_rank_file(gaston.results.file, "gaston_rank")

results <- join_results(r.results, spateo.results, gaston.results) %>%
  filter(
    .data$pattern == "peak",
    .data$peak.t.value == 1,
    .data$sd.r.value %in% c(0.5, 1, 2)
  )

if(nrow(results) == 0) {
  stop("No peak simulation rows found for peak.t = 1 and sd.r in 0.5, 1, 2.")
}

xy <- read_csv(file.path(data.dir, "granule_xy.csv"), show_col_types = FALSE)
t.coord <- read_csv(file.path(data.dir, "granule_t_coords.csv"), show_col_types = FALSE)[[1]]
r.coord <- read_csv(file.path(data.dir, "granule_r_coords.csv"), show_col_types = FALSE)[[1]]

rank.cols <- c(
  "mgam_rank",
  "spark_rank",
  "nnsvg_rank",
  "heartsvg_rank",
  "spateo_rank",
  "gaston_rank"
)
rank.cols <- intersect(rank.cols, colnames(results))

if(length(rank.cols) == 0) {
  stop("No rank columns were found in the available result files.")
}

method.labels <- c(
  mgam_rank = "MorphoGAM",
  spark_rank = "SPARK-X",
  nnsvg_rank = "nnSVG",
  heartsvg_rank = "HEARTSVG",
  spateo_rank = "Spateo",
  gaston_rank = "GASTON"
)

method.palette <- c(
  "MorphoGAM" = "#0072B2",
  "SPARK-X" = "#009E73",
  "nnSVG" = "#D55E00",
  "HEARTSVG" = "#CC79A7",
  "Spateo" = "#E69F00",
  "GASTON" = "#56B4E9"
)

plot.theme <- theme_bw(base_size = 11) +
  theme(
    axis.title = element_text(),
    axis.text = element_text(color = "black"),
    strip.background = element_rect(fill = "grey92", color = "grey70"),
    strip.text = element_text(),
    legend.position = "none",
    plot.title = element_text(hjust = 0),
    plot.subtitle = element_text(color = "grey30"),
    plot.margin = margin(4, 4, 4, 4)
  )

rescale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if(diff(rng) == 0) {
    return(rep(0, length(x)))
  }
  (x - rng[1]) / diff(rng)
}

peak_signal_from_simulator <- function(t, peak.t.value = 1) {
  T.basis <- ns(t, df = 10)

  set.seed(1)
  v <- rep(0, ncol(T.basis))
  v[sample(2:(ncol(T.basis) - 1), 1)] <- 1

  alpha1 <- rnorm(1, mean = -1)
  eta <- as.numeric(alpha1 + T.basis %*% (peak.t.value * v))
  rescale01(eta)
}

schematic.df <- tibble(
  xcoord = xy$xcoord,
  ycoord = xy$ycoord,
  gene1 = peak_signal_from_simulator(t.coord, peak.t.value = 1),
  confounder = rescale01(r.coord)
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

make_rank_panel <- function(results) {
  method.levels <- unname(method.labels[rank.cols])

  plot.df <- results %>%
    select(pattern, peak.t.value, sd.r.value, replicate, all_of(rank.cols)) %>%
    pivot_longer(
      cols = all_of(rank.cols),
      names_to = "method",
      values_to = "rank"
    ) %>%
    mutate(
      method = recode(.data$method, !!!method.labels),
      method = factor(.data$method, levels = method.levels),
      sd.r.value = factor(
        .data$sd.r.value,
        levels = c(0.5, 1, 2),
        labels = c("Low r", "Medium r", "High r")
      )
    ) %>%
    filter(!is.na(.data$rank), .data$rank > 0)

  if(nrow(plot.df) == 0) {
    stop("Rank columns were present, but all ranks were missing.")
  }

  ggplot(plot.df, aes(x = method, y = rank, fill = method)) +
    geom_boxplot(width = 0.62, outlier.shape = NA, linewidth = 0.35, alpha = 0.88) +
    geom_point(
      aes(color = method),
      position = position_jitter(width = 0.12, height = 0, seed = 1),
      size = 1.45,
      alpha = 0.85
    ) +
    facet_wrap(~ sd.r.value, nrow = 1) +
    scale_fill_manual(values = method.palette[method.levels]) +
    scale_color_manual(values = method.palette[method.levels]) +
    scale_y_sqrt(
      breaks = pretty_breaks(n = 5),
      expand = expansion(mult = c(0.05, 0.12))
    ) +
    labs(
      title = NULL,
      subtitle = NULL,
      x = NULL,
      y = "Rank of gene 1"
    ) +
    plot.theme +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1),
      panel.spacing.x = unit(0.9, "lines")
    )
}

p.gene1 <- make_spatial_panel(
  schematic.df,
  value_col = "gene1",
  title = "Gene 1: Localized signal",
  legend_title = "Signal"
)

p.confounder <- make_spatial_panel(
  schematic.df,
  value_col = "confounder",
  title = "Confounding genes 2 - G:\nlinear radial signal",
  legend_title = "Signal"
)

p.rank <- make_rank_panel(results)

p <- ggarrange(
  ggarrange(p.gene1, p.confounder, nrow = 1, ncol = 2),
  p.rank,
  nrow = 2,
  heights = c(1, 1.15)
)

pdf.file <- file.path(plot.dir, "granule_sim_results_peak.pdf")
png.file <- file.path(plot.dir, "granule_sim_results_peak.png")

ggsave(pdf.file, p, width = 9.11 * 1.1, height = 6, units = "in", bg = "white")
ggsave(png.file, p, width = 9.11 * 1.1, height = 6, units = "in", bg = "white")

message("Wrote ", normalizePath(pdf.file))
message("Wrote ", normalizePath(png.file))
