setwd(here::here("analysis/slideseq_ca3/code"))

library(here)
library(magick)

sim.plot.dir <- here::here("analysis", "ca3_newsim", "plots")

gradient.file <- file.path(sim.plot.dir, "ca3_sim_results_gradient.png")
nonparam.file <- file.path(sim.plot.dir, "ca3_sim_results_nonparam.png")

if(!file.exists(gradient.file)) {
  stop("Gradient simulation plot not found: ", gradient.file)
}

if(!file.exists(nonparam.file)) {
  stop("Non-parametric simulation plot not found: ", nonparam.file)
}

add_panel_label <- function(img, label) {
  image_annotate(
    img,
    text = label,
    gravity = "northwest",
    location = "+38+28",
    size = 82,
    weight = 700,
    color = "black"
  )
}

p.gradient <- image_read(gradient.file) |>
  add_panel_label("a")

p.nonparam <- image_read(nonparam.file) |>
  add_panel_label("b")

supp.plot <- image_append(c(p.gradient, p.nonparam), stack = TRUE)

out.file <- "../plots/ca3_sim_results_gradient_nonparam_supp.png"
image_write(supp.plot, path = out.file, format = "png")

message("Wrote ", normalizePath(out.file, mustWork = FALSE))
