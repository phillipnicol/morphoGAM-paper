#!/usr/bin/env Rscript

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

pkg <- "phillipnicol/MorphoGAM@rev_plan"

renv::install(pkg, rebuild = TRUE)
renv::snapshot(prompt = FALSE, force=TRUE)

desc <- utils::packageDescription("MorphoGAM")

dir.create("results/metadata",
           recursive = TRUE,
           showWarnings = FALSE)

meta <- c(
  paste0("package: ", desc$Package),
  paste0("remote: ",
         desc$RemoteUsername %||% "",
         "/",
         desc$RemoteRepo %||% ""),
  paste0("branch: ", desc$RemoteRef %||% ""),
  paste0("commit: ", desc$RemoteSha %||% ""),
  paste0("version: ", desc$Version),
  paste0("updated_at: ",
         format(Sys.time(), tz = "UTC"))
)

writeLines(
  meta,
  "results/metadata/morphogam-version.yml"
)