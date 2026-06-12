

#load("../data/ca3_curve.RData")
#out <- fit

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

data.dir <- file.path("..", "data")

cli.args <- commandArgs(trailingOnly = TRUE)

get_cli_value <- function(flag, default = NULL) {
  ix <- match(flag, cli.args)
  if(is.na(ix) || ix == length(cli.args)) {
    return(default)
  }
  cli.args[ix + 1]
}

has_cli_flag <- function(flag) {
  flag %in% cli.args
}

parse_numeric_vector <- function(x) {
  as.numeric(strsplit(x, ",", fixed = TRUE)[[1]])
}

parse_character_vector <- function(x) {
  strsplit(x, ",", fixed = TRUE)[[1]]
}

if(has_cli_flag("--help")) {
  cat(
    "Usage: Rscript run_rmodels.R [options]\n",
    "\n",
    "Options:\n",
    "  --replicates N        Number of replicates to run, starting at 1. Default: 10\n",
    "  --peak-t X,Y          Comma-separated peak.t values. Default: 1\n",
    "  --sd-r X,Y            Comma-separated sd.r values. Default: 0.5,1,2\n",
    "  --patterns A,B        Comma-separated patterns. Default: peak,gradient,nonparam\n",
    "  --output-file PATH    Results CSV path. Default: ../data/ca3_rmodel_results.csv\n",
    "  --dry-run             Print the simulation grid and exit.\n",
    sep = ""
  )
  quit(save = "no", status = 0)
}

M <- as.integer(get_cli_value("--replicates", "10"))
peak.t <- parse_numeric_vector(get_cli_value("--peak-t", "1"))
sd.r <- parse_numeric_vector(get_cli_value("--sd-r", "0.5,1,2"))
patterns <- parse_character_vector(get_cli_value("--patterns", "peak,gradient,nonparam"))
results.file <- get_cli_value("--output-file", file.path(data.dir, "ca3_rmodel_results.csv"))

params <- expand.grid(
  pattern = patterns,
  peak.t.value = peak.t,
  sd.r.value = sd.r,
  replicate = seq_len(M),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

if(has_cli_flag("--dry-run")) {
  print(params)
  message("Would process ", nrow(params), " datasets.")
  quit(save = "no", status = 0)
}


library(MorphoGAM)



spe <- STexampleData::SlideSeqV2_mouseHPC()

ixs <- which(spe$celltype == "CA3") #subset to CA3

xy <- spatialCoords(spe)[ixs,]
Y <- counts(spe)[,ixs]

xy.dist <- as.matrix(dist(xy))
knn <- 20
prune.outlier <- 3
nnk <- apply(xy.dist, 1, function(x) sort(x)[knn+1])
outlier <- which(nnk > (prune.outlier)*median(nnk))

xy <- xy[-outlier,]; Y <- Y[,-outlier]
df <- xy

fit <- CurveFinder(xy) #Fit the curve
out <- fit


set.seed(1)

format_label <- function(x) {
  gsub("\\.", "p", format(x, trim = TRUE, scientific = FALSE))
}

rank_gene <- function(gene_ids, gene = "1") {
  rank <- match(gene, gene_ids)
  if(is.na(rank)) NA_integer_ else rank
}

read_simulated_counts <- function(in.file) {
  if(!file.exists(in.file)) {
    stop("Missing simulated dataset: ", in.file)
  }

  Y_sim <- read.csv(in.file, row.names = 1, check.names = FALSE)
  Y_sim <- as.matrix(Y_sim)
  storage.mode(Y_sim) <- "numeric"
  rownames(Y_sim) <- as.character(seq_len(nrow(Y_sim)))
  Y_sim
}

run_morphogam <- function(Y_sim) {
  mgam <- MorphoGAM::MorphoGAM(
    Y = Y_sim,
    curve.fit = fit,
    design = y ~ s(t, bs = "cr", k = 20) + s(r),
    return.fx = TRUE
  )

  mgam.results <- mgam$results
  ord <- order(mgam.results$peak.t, decreasing = TRUE, na.last = TRUE)
  gene.ids <- rownames(mgam.results)[ord]

  list(
    rank = rank_gene(gene.ids),
    score = mgam.results["1", "peak.t"]
  )
}

run_sparkx <- function(Y_sim) {
  spark.res <- SPARK::sparkx(
    count_in = Y_sim,
    locus_in = xy
  )

  spark.results <- spark.res$res_mtest
  ord <- order(spark.results$adjustedPval, decreasing = FALSE, na.last = TRUE)
  gene.ids <- rownames(spark.results)[ord]

  list(
    rank = rank_gene(gene.ids),
    score = spark.results["1", "adjustedPval"]
  )
}

run_nnsvg <- function(Y_sim) {
  library(nnSVG)

  lib.size <- colSums(Y_sim)
  lib.size[lib.size == 0] <- 1
  logCPM <- log(sweep(Y_sim, MARGIN = 2, STATS = lib.size, FUN = "/") + 1)

  nn.svg <- nnSVG(
    input = logCPM,
    spatial_coords = xy,
    verbose = TRUE
  )

  if(nrow(nn.svg) < 1 || ncol(nn.svg) < 12) {
    stop("nnSVG result does not have entry [1, 12].")
  }

  gene1.rank <- nn.svg[1, 12]
  score.col <- colnames(nn.svg)[12]

  list(rank = gene1.rank, score = gene1.rank, score_col = score.col)
}

run_heartsvg <- function(Y_sim) {
  hsvg.mat <- data.frame(
    row = xy[, 1],
    col = xy[, 2],
    t(Y_sim),
    check.names = FALSE
  )
  colnames(hsvg.mat) <- c("row", "col", paste0("gene", seq_len(nrow(Y_sim))))

  heartsvg.results <- HEARTSVG::heartsvg(hsvg.mat, scale = TRUE)
  gene.ids <- as.character(heartsvg.results$gene)
  rank <- rank_gene(gene.ids, gene = "gene1")

  score.cols <- setdiff(colnames(heartsvg.results), "gene")
  numeric.cols <- score.cols[vapply(heartsvg.results[score.cols], is.numeric, logical(1))]
  score.col <- if(length(numeric.cols) > 0) numeric.cols[1] else NA_character_
  score <- if(!is.na(rank) && !is.na(score.col)) {
    heartsvg.results[rank, score.col]
  } else {
    NA_real_
  }

  list(
    rank = rank,
    score = score,
    score_col = score.col
  )
}

empty_result <- function(pattern, peak.t.value, sd.r.value, replicate, in.file) {
  data.frame(
    pattern = pattern,
    peak.t.value = peak.t.value,
    sd.r.value = sd.r.value,
    replicate = replicate,
    input_file = in.file,
    mgam_rank = NA_integer_,
    mgam_peak_t = NA_real_,
    mgam_error = NA_character_,
    spark_rank = NA_integer_,
    spark_adjusted_pval = NA_real_,
    spark_error = NA_character_,
    nnsvg_rank = NA_integer_,
    nnsvg_score = NA_real_,
    nnsvg_score_col = NA_character_,
    nnsvg_error = NA_character_,
    heartsvg_rank = NA_integer_,
    heartsvg_score = NA_real_,
    heartsvg_score_col = NA_character_,
    heartsvg_error = NA_character_,
    stringsAsFactors = FALSE
  )
}

run_one_dataset <- function(pattern, peak.t.value, sd.r.value, replicate) {
  peak.t.label <- format_label(peak.t.value)
  sd.r.label <- format_label(sd.r.value)
  pattern.dir <- file.path(data.dir, paste0("ca3_simulation_", pattern))
  in.file <- file.path(
    pattern.dir,
    sprintf(
      "ca3_simulated_dataset_peakt_%s_sd.r_%s_rep_%02d.csv",
      peak.t.label,
      sd.r.label,
      replicate
    )
  )

  message(
    "Testing pattern = ", pattern,
    ", peak.t = ", peak.t.value,
    ", sd.r = ", sd.r.value,
    ", replicate = ", replicate
  )

  result <- empty_result(pattern, peak.t.value, sd.r.value, replicate, in.file)
  Y_sim <- read_simulated_counts(in.file)

  mgam.out <- tryCatch(run_morphogam(Y_sim), error = function(e) e)
  if(inherits(mgam.out, "error")) {
    result$mgam_error <- conditionMessage(mgam.out)
  } else {
    result$mgam_rank <- mgam.out$rank
    result$mgam_peak_t <- mgam.out$score
  }

  spark.out <- tryCatch(run_sparkx(Y_sim), error = function(e) e)
  if(inherits(spark.out, "error")) {
    result$spark_error <- conditionMessage(spark.out)
  } else {
    result$spark_rank <- spark.out$rank
    result$spark_adjusted_pval <- spark.out$score
  }

  nnsvg.out <- tryCatch(run_nnsvg(Y_sim), error = function(e) e)
  if(inherits(nnsvg.out, "error")) {
    result$nnsvg_error <- conditionMessage(nnsvg.out)
  } else {
    result$nnsvg_rank <- nnsvg.out$rank
    result$nnsvg_score <- nnsvg.out$score
    result$nnsvg_score_col <- nnsvg.out$score_col
  }

  heartsvg.out <- tryCatch(run_heartsvg(Y_sim), error = function(e) e)
  if(inherits(heartsvg.out, "error")) {
    result$heartsvg_error <- conditionMessage(heartsvg.out)
  } else {
    result$heartsvg_rank <- heartsvg.out$rank
    result$heartsvg_score <- heartsvg.out$score
    result$heartsvg_score_col <- heartsvg.out$score_col
  }

  result
}

results <- vector("list", nrow(params))

for(i in seq_len(nrow(params))) {
  results[[i]] <- run_one_dataset(
    pattern = params$pattern[i],
    peak.t.value = params$peak.t.value[i],
    sd.r.value = params$sd.r.value[i],
    replicate = params$replicate[i]
  )

  results.df <- do.call(rbind, results[seq_len(i)])
  write.csv(results.df, results.file, row.names = FALSE)
  message("Checkpointed ", i, " / ", nrow(params), " rows to ", results.file)
}

message("Wrote final R model results to ", normalizePath(results.file))
