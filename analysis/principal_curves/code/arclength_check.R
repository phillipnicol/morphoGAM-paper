# swissroll_curvefinder_geodesic_mse_finite.R
#
# Purpose:
#   Compare CurveFinder-style kNN graph geodesic distances to true
#   arclength distances on a known swiss-roll-like curve.
#
# Output:
#   ../results/swissroll_curvefinder_geodesic_mse_finite_raw.csv
#   ../results/swissroll_curvefinder_geodesic_mse_finite_summary.csv
#   ../plots/swissroll_curvefinder_geodesic_mse_finite.png

library(tidyverse)
library(igraph)
library(dimRed)
library(viridis)

set.seed(1)

# -----------------------------
# User settings
# -----------------------------

n_grid <- c(50, 100, 200, 500, 1000)
k_grid <- 3:10
n_reps <- 10
noise_sd <- 0.01

out_dir <- "../results"
plot_dir <- "../plots"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# True swiss-roll curve
# -----------------------------
#
# Curve:
#   f(t) = r(t) (cos(10t), sin(10t))
#   r(t) = 0.5 + 0.5t
#
# Speed:
#   ||f'(t)|| = sqrt( r'(t)^2 + (r(t) theta'(t))^2 )
#             = sqrt(0.5^2 + (10 r(t))^2)

speed_fun <- function(t) {
  r <- 0.5 + 0.5 * t
  sqrt(0.5^2 + (10 * r)^2)
}

make_arclength_interpolator <- function(grid_size = 50000) {
  grid <- seq(0, 1, length.out = grid_size)
  speed <- speed_fun(grid)

  arc_grid <- c(
    0,
    cumsum(diff(grid) * (head(speed, -1) + tail(speed, -1)) / 2)
  )

  list(
    grid = grid,
    arc_grid = arc_grid,
    total_length = max(arc_grid)
  )
}

arc_obj <- make_arclength_interpolator()

true_arclength <- function(t) {
  approx(
    x = arc_obj$grid,
    y = arc_obj$arc_grid,
    xout = t,
    rule = 2
  )$y
}

simulate_swissroll <- function(n, noise_sd = 0.01) {
  t <- runif(n = n, min = 0, max = 1)
  r <- 0.5 + 0.5 * t

  xy_noiseless <- cbind(
    x = r * cos(10 * t),
    y = r * sin(10 * t)
  )

  xy <- xy_noiseless + matrix(
    rnorm(2 * n, mean = 0, sd = noise_sd),
    nrow = n,
    ncol = 2
  )

  list(
    xy = xy,
    t = t,
    s_true = true_arclength(t)
  )
}

# -----------------------------
# CurveFinder-style geodesic
# -----------------------------
#
# This mirrors the relevant calculation inside MorphoGAM::CurveFinder():
#
#   knng <- dimRed:::makeKNNgraph(x = xy, k = knn, eps = 0)
#   dG <- igraph::distances(knng, algorithm = "dijkstra")

compute_curvefinder_geodesic <- function(xy, k) {
  knng <- dimRed:::makeKNNgraph(
    x = xy,
    k = k,
    eps = 0
  )

  igraph::distances(
    knng,
    algorithm = "dijkstra"
  )
}

# -----------------------------
# Evaluate one simulation
# -----------------------------

evaluate_one <- function(n, k, rep_id, noise_sd) {
  dat <- simulate_swissroll(n = n, noise_sd = noise_sd)

  dG <- compute_curvefinder_geodesic(
    xy = dat$xy,
    k = k
  )

  # True pairwise arclength distance matrix.
  d_true <- abs(outer(dat$s_true, dat$s_true, "-"))

  # Normalize by total true arclength so MSE is unitless.
  dG_norm <- dG / arc_obj$total_length
  d_true_norm <- d_true / arc_obj$total_length

  # Use upper triangle to avoid duplicate pairs and diagonal zeros.
  upper <- upper.tri(d_true_norm)

  dG_vec <- dG_norm[upper]
  d_true_vec <- d_true_norm[upper]

  finite <- is.finite(dG_vec)

  mse_finite <- mean(
    (dG_vec[finite] - d_true_vec[finite])^2
  )

  tibble(
    n = n,
    k = k,
    rep = rep_id,
    mse_finite = mse_finite,
    prop_finite_pairs = mean(finite)
  )
}

# -----------------------------
# Run grid
# -----------------------------

grid <- tidyr::crossing(
  n = n_grid,
  k = k_grid,
  rep = seq_len(n_reps)
) |>
  dplyr::filter(k < n)

res <- purrr::pmap_dfr(
  grid,
  function(n, k, rep) {
    message("Running n = ", n, ", k = ", k, ", rep = ", rep)

    evaluate_one(
      n = n,
      k = k,
      rep_id = rep,
      noise_sd = noise_sd
    )
  }
)

summary_res <- res |>
  dplyr::group_by(n, k) |>
  dplyr::summarize(
    mean_mse_finite = mean(mse_finite, na.rm = TRUE),
    se_mse_finite = sd(mse_finite, na.rm = TRUE) /
      sqrt(sum(!is.na(mse_finite))),
    mean_prop_finite_pairs = mean(prop_finite_pairs, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(
  res,
  file.path(out_dir, "swissroll_curvefinder_geodesic_mse_finite_raw.csv")
)

readr::write_csv(
  summary_res,
  file.path(out_dir, "swissroll_curvefinder_geodesic_mse_finite_summary.csv")
)

# -----------------------------
# Plot
# -----------------------------

label_cutoff <- stats::quantile(
  summary_res$mean_mse_finite,
  probs = 0.75,
  na.rm = TRUE
)

p_mse_finite <- summary_res |>
  dplyr::mutate(
    label_color = ifelse(mean_mse_finite >= label_cutoff, "white", "grey15")
  ) |>
  ggplot(aes(x = factor(k), y = factor(n), fill = mean_mse_finite)) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(
    aes(label = sprintf("%.3f", mean_mse_finite), color = label_color),
    size = 3.2
  ) +
  scale_color_identity() +
  scale_fill_gradient(
    low = "grey92",
    high = "firebrick3",
    trans = "sqrt",
    name = "MSE"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10, color = "grey20"),
    plot.title = element_text(size = 13, hjust = 0),
    plot.subtitle = element_text(size = 10, color = "grey35", hjust = 0),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.key.height = unit(0.7, "cm"),
    plot.margin = margin(8, 10, 8, 8)
  ) +
  labs(
    x = "Number of nearest neighbors (k)",
    y = "Number of sampled points (n)",
    title = "MSE between graph shortest paths and true arclengths"
  )

ggsave(
  filename = file.path(plot_dir, "swissroll_curvefinder_geodesic_mse_finite.png"),
  plot = p_mse_finite,
  width = 7,
  height = 4.5,
  units = "in",
  dpi = 300
)