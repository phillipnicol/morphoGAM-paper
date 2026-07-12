setwd(here::here("analysis/merfish_swissroll/code"))

library(tidyverse)
library(ggpubr)
library(Matrix)

dir.create("../plots", showWarnings = FALSE, recursive = TRUE)

load_fit <- function(path) {
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  env$fit
}

plot_curve <- function(fit, title) {
  curve_df <- NULL
  if (!is.null(fit$curve.plot$layers[[2]]$data) &&
      is.data.frame(fit$curve.plot$layers[[2]]$data)) {
    curve_df <- fit$curve.plot$layers[[2]]$data
  }

  if (is.null(curve_df)) {
    curve_df <- fit$xyt |>
      arrange(t) |>
      transmute(x = f1, y = f2, color = t)
  }

  ggplot(fit$xyt, aes(x = x, y = y)) +
    geom_point(color = "grey", alpha = 0.5, size = 1.5) +
    geom_path(
      data = curve_df,
      aes(x = x, y = y, color = color),
      linewidth = 1,
      inherit.aes = FALSE
    ) +
    scale_color_gradient(low = "navyblue",
                         high = "firebrick1",
                         limits = c(0, 1)) +
    theme_bw() +
    labs(title = title, color = "t")
}

plot_counts_for_roll <- function(Y, fit, res, intervals, gene.names) {
  stopifnot(ncol(Y) == nrow(fit$xyt))

  nmed <- median(Matrix::colSums(Y))
  depth.scale <- Matrix::colSums(Y) / nmed
  indic <- Reduce(
    `|`,
    lapply(intervals, function(x) fit$xyt$t > x[1] & fit$xyt$t < x[2])
  )

  observed <- sweep(as.matrix(Y[gene.names, , drop = FALSE]),
                    MARGIN = 2,
                    STATS = depth.scale,
                    FUN = "/")

  fitted <- matrix(NA_real_,
                   nrow = length(gene.names),
                   ncol = ncol(Y),
                   dimnames = list(gene.names, colnames(Y)))
  fitted[, !indic] <- nmed * exp(res[gene.names, 2])
  fitted[, indic] <- nmed * exp(res[gene.names, 2] + res[gene.names, 1])

  df <- t(observed) |>
    as.data.frame() |>
    mutate(t = fit$xyt$t) |>
    pivot_longer(-t)

  df.fit <- t(fitted) |>
    as.data.frame() |>
    mutate(t = fit$xyt$t) |>
    pivot_longer(-t)

  df$fitted <- df.fit$value
  df
}

plot_ulcerated_regions <- function(fit, intervals) {
  indic <- Reduce(
    `|`,
    lapply(intervals, function(x) fit$xyt$t > x[1] & fit$xyt$t < x[2])
  )

  df <- data.frame(
    x = fit$xyt$x,
    y = fit$xyt$y,
    color = ifelse(indic, "Ulcerated", "Non-ulcerated")
  )

  ggplot(df, aes(x = x, y = y, color = color)) +
    geom_point(size = 0.25) +
    theme_bw() +
    scale_color_manual(values = c("Ulcerated" = "orange",
                                  "Non-ulcerated" = "grey90")) +
    labs(color = "") +
    xlab("x") +
    ylab("y")
}

count_plot <- function(df, intervals) {
  p <- ggplot(df, aes(x = t, y = value)) +
    theme_bw() +
    geom_point(size = 0.25) +
    ylab("Depth-normalized count") +
    geom_line(aes(x = t, y = fitted, color = name)) +
    facet_wrap(~name, ncol = 1, scales = "free_y") +
    labs(color = "Gene") +
    scale_y_sqrt() +
    guides(color = "none")

  for (interval in intervals) {
    p <- p +
      annotate("rect",
               xmin = interval[1],
               xmax = interval[2],
               ymin = 0,
               ymax = Inf,
               fill = "lightgrey",
               alpha = 0.5)
  }

  p
}

### PLOT STARTS HERE

res1 <- readRDS("../data/roll1_intercept_model.RDS")
res2 <- readRDS("../data/roll2_intercept_model.RDS")
res3 <- readRDS("../data/roll3_intercept_model.RDS")

beta.combine <- data.frame(
  gene = rownames(res1),
  beta1 = res1[, 1],
  beta2 = res2[, 1],
  beta3 = res3[, 1],
  pval = pmax(res1[, 3], res2[, 3], res3[, 3])
)
beta.combine$pval <- p.adjust(beta.combine$pval, method = "BH")

gene.names <- c("Mmp3", "Timp1", "Mmp10", "Il11", "Grem1")
gene_names <- readRDS("../data/gene_names.RDS")

roll1.intervals <- list(c(0.95, 0.97))
roll2.intervals <- list(c(0.25, 0.6), c(0.78, 0.8), c(0.88, 0.9))
roll3.intervals <- list(c(0.3, 0.6), c(0.63, 0.75), c(0.9, 0.94))

# 061923_D9_m2_Swiss is roll 1.
fit1 <- load_fit("../data/061923_D9_m2_Swiss_fibroblast.RData")
Y1 <- Matrix::readMM("../data/061923_D9_m2_Swiss_fibroblast_counts.mtx")
rownames(Y1) <- gene_names
df1 <- plot_counts_for_roll(Y1, fit1, res1, roll1.intervals, gene.names)
p1 <- count_plot(df1, roll1.intervals)
rm(Y1, df1)
gc()

# curve_d9_m5_080823 is roll 2.
fit2 <- load_fit("../data/curve_d9_m5_080823.RData")
Y2 <- Matrix::t(Matrix::readMM("../data/080823_d9_m5.mtx"))
rownames(Y2) <- gene_names
df2 <- plot_counts_for_roll(Y2, fit2, res2, roll2.intervals, gene.names)
p2 <- count_plot(df2, roll2.intervals)
rm(Y2, df2)
gc()

# 080823_D9_m13_Swiss_fibroblast is roll 3.
fit3 <- load_fit("../data/080823_D9_m13_Swiss_fibroblast.RData")
Y3 <- Matrix::readMM("../data/080823_D9_m13_Swiss_fibroblast_counts.mtx")
rownames(Y3) <- gene_names
df3 <- plot_counts_for_roll(Y3, fit3, res3, roll3.intervals, gene.names)
p3 <- count_plot(df3, roll3.intervals)
rm(Y3, df3)
gc()

p.swiss1 <- plot_curve(fit1, "Swiss roll 1") + guides(color = "none")
p.swiss2 <- plot_curve(fit2, "Swiss roll 2") + guides(color = "none")
p.swiss3 <- plot_curve(fit3, "Swiss roll 3") + guides(color = "none")

p.swiss.ulc1 <- plot_ulcerated_regions(fit1, roll1.intervals)
p.swiss.ulc2 <- plot_ulcerated_regions(fit2, roll2.intervals)
p.swiss.ulc3 <- plot_ulcerated_regions(fit3, roll3.intervals)

p <- ggarrange(
  ggarrange(p.swiss1, p.swiss2, p.swiss3, nrow = 1),
  ggarrange(p.swiss.ulc1, p.swiss.ulc2, p.swiss.ulc3,
            nrow = 1, common.legend = TRUE, legend = "top"),
  ggarrange(p1, p2, p3, nrow = 1),
  nrow = 3,
  heights = c(1, 1, 1.5),
  labels = c("a", "b", "c")
)

crush <- 0.5
ggsave(p, filename = "../plots/swissrolls_ulceratedgenes.png",
       height = crush * 6171, width = crush * 4911, units = "px",
       bg = "white")

# Just ulcerated.
p <- ggarrange(
  ggarrange(p.swiss1, p.swiss2, p.swiss3, nrow = 1),
  ggarrange(p.swiss.ulc1, p.swiss.ulc2, p.swiss.ulc3,
            nrow = 1, common.legend = TRUE, legend = "top"),
  nrow = 2
)

ggsave(p, filename = "../plots/ulcerated_swissrolls.png", bg = "white")

p <- ggarrange(p1, p2, p3, nrow = 1)

ggsave(p, filename = "../plots/ulcerated_iaf_markers.png", bg = "white")
