setwd(here::here("analysis/ca3_newsim/code"))

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
library(ggplot2)
library(patchwork)
library(readr)
library(scales)
library(tidyr)

data.dir <- file.path("..", "data")
plot.dir <- file.path("..", "plots")
dir.create(plot.dir, showWarnings = FALSE, recursive = TRUE)

full.results.file <- file.path(data.dir, "ca3_rmodel_results.csv")
smoke.results.file <- file.path(data.dir, "ca3_rmodel_results_smoke.csv")
results.file <- if(file.exists(full.results.file)) full.results.file else smoke.results.file

if(!file.exists(results.file)) {
  stop("No results file found. Expected ", full.results.file, " or ", smoke.results.file)
}

message("Reading simulation results from ", results.file)
results <- read_csv(results.file, show_col_types = FALSE)

spateo.results.file <- file.path(data.dir, "ca3_spateo_model_results.csv")
if(file.exists(spateo.results.file)) {
  message("Adding Spateo results from ", spateo.results.file)
  spateo.results <- read_csv(spateo.results.file, show_col_types = FALSE) %>%
    select(pattern, peak.t.value, sd.r.value, replicate, spateo_rank)

  results <- results %>%
    left_join(
      spateo.results,
      by = c("pattern", "peak.t.value", "sd.r.value", "replicate")
    )
} else {
  message("No Spateo results file found at ", spateo.results.file)
}

gaston.results.file <- file.path(data.dir, "ca3_gaston_model_results.csv")
if(file.exists(gaston.results.file)) {
  message("Adding GASTON results from ", gaston.results.file)
  gaston.results <- read_csv(gaston.results.file, show_col_types = FALSE) %>%
    select(pattern, peak.t.value, sd.r.value, replicate, gaston_rank)

  results <- results %>%
    left_join(
      gaston.results,
      by = c("pattern", "peak.t.value", "sd.r.value", "replicate")
    )
} else {
  message("No GASTON results file found at ", gaston.results.file)
}

xy <- read_csv(file.path(data.dir, "ca3_xy.csv"), show_col_types = FALSE)
t.coord <- read_csv(file.path(data.dir, "ca3_t_coords.csv"), show_col_types = FALSE)[[1]]
r.coord <- read_csv(file.path(data.dir, "ca3_r_coords.csv"), show_col_types = FALSE)[[1]]

rank.cols <- c(
  "mgam_rank",
  "spark_rank",
  "nnsvg_rank",
  "heartsvg_rank",
  "spateo_rank",
  "gaston_rank"
)

missing.rank.cols <- setdiff(rank.cols, colnames(results))
if(length(missing.rank.cols) > 0) {
  stop("Missing rank columns in ", results.file, ": ", paste(missing.rank.cols, collapse = ", "))
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

plot.theme <- theme_bw(base_size = 12) +
  theme(
    axis.title = element_text(),
    axis.text = element_text(color = "black"),
    strip.background = element_rect(fill = "grey92", color = "grey70"),
    strip.text = element_text(),
    legend.position = "none",
    plot.title = element_text(hjust = 0),
    plot.subtitle = element_text(color = "grey30"),
    plot.margin = margin(6, 6, 6, 6)
  )

rescale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if(diff(rng) == 0) {
    return(rep(0, length(x)))
  }
  (x - rng[1]) / diff(rng)
}

gene1_signal <- function(pattern, t) {
  if(pattern == "peak") {
    signal <- exp(-0.5 * ((t - 0.5) / 0.12)^2)
  } else if(pattern == "gradient") {
    signal <- t
  } else if(pattern == "nonparam") {
    signal <- 0.55 * sin(2 * pi * t) + 0.35 * cos(4 * pi * t) + 0.25 * t
  } else {
    stop("Unknown pattern: ", pattern)
  }

  rescale01(signal)
}

pattern_title <- function(pattern) {
  dplyr::case_when(
    pattern == "peak" ~ "Localized",
    pattern == "gradient" ~ "Gradient",
    pattern == "nonparam" ~ "Non-parametric",
    TRUE ~ pattern
  )
}

schematic_df <- function(pattern) {
  tibble(
    xcoord = xy$xcoord,
    ycoord = xy$ycoord,
    gene1 = gene1_signal(pattern, t.coord),
    confounder = rescale01(r.coord)
  )
}

make_spatial_panel <- function(df, value_col, title, legend_title) {
  ggplot(df, aes(x = xcoord, y = ycoord, color = .data[[value_col]])) +
    geom_point(size = 0.65, alpha = 0.95) +
    coord_equal() +
    scale_color_viridis_c(option = "magma", name = legend_title) +
    labs(title = title, x = "x", y = "y") +
    guides(color = guide_colorbar(barheight = unit(32, "pt"))) +
    theme_bw(base_size = 12) +
    theme(
      axis.title = element_text(),
      axis.text = element_text(color = "black"),
      legend.position = "right",
      legend.title = element_text(),
      plot.title = element_text(hjust = 0),
      plot.margin = margin(6, 6, 6, 6)
    )
}

make_rank_panel <- function(results, pattern) {
  plot.df <- results %>%
    filter(.data$pattern == .env$pattern) %>%
    select(pattern, peak.t.value, sd.r.value, replicate, all_of(rank.cols)) %>%
    pivot_longer(
      cols = all_of(rank.cols),
      names_to = "method",
      values_to = "rank"
    ) %>%
    mutate(
      method = recode(.data$method, !!!method.labels),
      method = factor(.data$method, levels = unname(method.labels)),
      sd.r.value = factor(
        .data$sd.r.value,
        levels = c(0.5, 1, 2),
        labels = c("Low r", "Medium r", "High r")
      )
    ) %>%
    filter(!is.na(.data$rank), .data$rank > 0)

  ggplot(plot.df, aes(x = method, y = rank, fill = method)) +
    geom_boxplot(width = 0.62, outlier.shape = NA, linewidth = 0.35, alpha = 0.88) +
    geom_point(
      aes(color = method),
      position = position_jitter(width = 0.12, height = 0, seed = 1),
      size = 1.8,
      alpha = 0.85
    ) +
    facet_wrap(~ sd.r.value, nrow = 1) +
    scale_fill_manual(values = method.palette) +
    scale_color_manual(values = method.palette) +
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
      panel.spacing.x = unit(1.1, "lines")
    )
}

make_pattern_plot <- function(pattern) {
  schem <- schematic_df(pattern)

  p.gene1 <- make_spatial_panel(
    schem,
    value_col = "gene1",
    title = paste0(pattern_title(pattern), " signal"),
    legend_title = "Signal"
  )

  p.confounder <- make_spatial_panel(
    schem,
    value_col = "confounder",
    title = "Confounding genes 2 - G: linear radial signal",
    legend_title = "Signal"
  )

  p.rank <- make_rank_panel(results, pattern)

  ((p.gene1 | p.confounder) / p.rank) +
    plot_layout(heights = c(1, 1.15))
}

patterns <- c("peak", "gradient", "nonparam")

for(pattern in patterns) {
  p <- make_pattern_plot(pattern)

  pdf.file <- file.path(plot.dir, paste0("ca3_sim_results_", pattern, ".pdf"))
  png.file <- file.path(plot.dir, paste0("ca3_sim_results_", pattern, ".png"))

  ggsave(pdf.file, p, width = 9.5, height = 7.2, units = "in", device = cairo_pdf)
  ggsave(png.file, p, width = 9.5, height = 7.2, units = "in", dpi = 300)

  message("Wrote ", pdf.file)
  message("Wrote ", png.file)
}
