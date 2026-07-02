

#load("../data/granule_curve.RData")
#out <- fit

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
    "  --methods A,B         Comma-separated methods to run: mgam,spark,nnsvg,heartsvg. Default: all\n",
    "  --output-file PATH    Results CSV path. Default: ../data/granule_rmodel_results.csv\n",
    "  --dry-run             Print the simulation grid and exit.\n",
    sep = ""
  )
  quit(save = "no", status = 0)
}

M <- as.integer(get_cli_value("--replicates", "10"))
peak.t <- parse_numeric_vector(get_cli_value("--peak-t", "1"))
sd.r <- parse_numeric_vector(get_cli_value("--sd-r", "0.5,1,2"))
patterns <- parse_character_vector(get_cli_value("--patterns", "peak,gradient,nonparam"))
methods <- parse_character_vector(get_cli_value("--methods", "mgam,spark,nnsvg,heartsvg"))
results.file <- get_cli_value("--output-file", file.path(data.dir, "granule_rmodel_results.csv"))
if(grepl("^analysis/", results.file)) {
  results.file <- here::here(results.file)
}
dir.create(dirname(results.file), showWarnings = FALSE, recursive = TRUE)

valid.methods <- c("mgam", "spark", "nnsvg", "heartsvg")
unknown.methods <- setdiff(methods, valid.methods)
if(length(unknown.methods) > 0) {
  stop("Unknown --methods value(s): ", paste(unknown.methods, collapse = ", "))
}

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
  message("Methods: ", paste(methods, collapse = ", "))
  message("Would process ", nrow(params), " datasets.")
  quit(save = "no", status = 0)
}

message("Methods selected: ", paste(methods, collapse = ", "))
message("Output file: ", results.file)

library(MorphoGAM)


prepare_granule_curve <- function() {
  xy <- readRDS(here::here("analysis", "data", "granule_slideseq.RDS"))
  xy <- as.matrix(xy)

  xy.dist <- as.matrix(dist(xy))
  nnk <- apply(xy.dist, 1, function(x) sort(x)[6])
  outlier <- which(nnk > 2 * median(nnk))
  if(length(outlier) > 0) {
    xy <- xy[-outlier, , drop = FALSE]
  }

  colnames(xy) <- c("xcoord", "ycoord")
  #Load manual curve with t and r coord
  load("../../principal_curves/data/granule_ground_truth.Rda")
  fit <- ground_truth

  list(xy = xy, fit = fit)
}

granule <- prepare_granule_curve()
xy <- granule$xy
fit <- granule$fit

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
    design = y ~ s(t) + s(r),
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
  heartsvg.xy <- apply(xy, 2, function(z) {
    rng <- range(z, na.rm = TRUE)
    if(diff(rng) == 0) {
      return(rep(0, length(z)))
    }
    (z - rng[1]) / diff(rng) * 1000
  })

  hsvg.mat <- data.frame(
    row = heartsvg.xy[, 1],
    col = heartsvg.xy[, 2],
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

input_file_for <- function(pattern, peak.t.value, sd.r.value, replicate) {
  peak.t.label <- format_label(peak.t.value)
  sd.r.label <- format_label(sd.r.value)
  pattern.dir <- file.path(data.dir, paste0("granule_simulation_", pattern))
  file.path(
    pattern.dir,
    sprintf(
      "granule_simulated_dataset_peakt_%s_sd.r_%s_rep_%02d.csv",
      peak.t.label,
      sd.r.label,
      replicate
    )
  )
}

existing.results <- if(file.exists(results.file)) {
  read.csv(results.file, check.names = FALSE, stringsAsFactors = FALSE)
} else {
  data.frame()
}

existing_result <- function(pattern, peak.t.value, sd.r.value, replicate, in.file) {
  result <- empty_result(pattern, peak.t.value, sd.r.value, replicate, in.file)

  if(nrow(existing.results) == 0) {
    return(result)
  }

  ix <- which(
    existing.results$pattern == pattern &
      as.numeric(existing.results$peak.t.value) == peak.t.value &
      as.numeric(existing.results$sd.r.value) == sd.r.value &
      as.integer(existing.results$replicate) == replicate
  )

  if(length(ix) == 0) {
    return(result)
  }

  old.result <- existing.results[ix[1], , drop = FALSE]
  common.cols <- intersect(colnames(result), colnames(old.result))
  result[common.cols] <- old.result[common.cols]
  result$input_file <- in.file
  result
}

run_one_dataset <- function(pattern, peak.t.value, sd.r.value, replicate) {
  in.file <- input_file_for(pattern, peak.t.value, sd.r.value, replicate)

  message(
    "Testing pattern = ", pattern,
    ", peak.t = ", peak.t.value,
    ", sd.r = ", sd.r.value,
    ", replicate = ", replicate
  )

  result <- existing_result(pattern, peak.t.value, sd.r.value, replicate, in.file)
  Y_sim <- read_simulated_counts(in.file)

  if("mgam" %in% methods) {
    result$mgam_rank <- NA_integer_
    result$mgam_peak_t <- NA_real_
    result$mgam_error <- NA_character_

    mgam.out <- tryCatch(run_morphogam(Y_sim), error = function(e) e)
    if(inherits(mgam.out, "error")) {
      result$mgam_error <- conditionMessage(mgam.out)
    } else {
      result$mgam_rank <- mgam.out$rank
      result$mgam_peak_t <- mgam.out$score
    }
  }

  if("spark" %in% methods) {
    result$spark_rank <- NA_integer_
    result$spark_adjusted_pval <- NA_real_
    result$spark_error <- NA_character_

    spark.out <- tryCatch(run_sparkx(Y_sim), error = function(e) e)
    if(inherits(spark.out, "error")) {
      result$spark_error <- conditionMessage(spark.out)
    } else {
      result$spark_rank <- spark.out$rank
      result$spark_adjusted_pval <- spark.out$score
    }
  }

  if("nnsvg" %in% methods) {
    result$nnsvg_rank <- NA_integer_
    result$nnsvg_score <- NA_real_
    result$nnsvg_score_col <- NA_character_
    result$nnsvg_error <- NA_character_

    nnsvg.out <- tryCatch(run_nnsvg(Y_sim), error = function(e) e)
    if(inherits(nnsvg.out, "error")) {
      result$nnsvg_error <- conditionMessage(nnsvg.out)
    } else {
      result$nnsvg_rank <- nnsvg.out$rank
      result$nnsvg_score <- nnsvg.out$score
      result$nnsvg_score_col <- nnsvg.out$score_col
    }
  }

  if("heartsvg" %in% methods) {
    result$heartsvg_rank <- NA_integer_
    result$heartsvg_score <- NA_real_
    result$heartsvg_score_col <- NA_character_
    result$heartsvg_error <- NA_character_

    heartsvg.out <- tryCatch(run_heartsvg(Y_sim), error = function(e) e)
    if(inherits(heartsvg.out, "error")) {
      result$heartsvg_error <- conditionMessage(heartsvg.out)
    } else {
      result$heartsvg_rank <- heartsvg.out$rank
      result$heartsvg_score <- heartsvg.out$score
      result$heartsvg_score_col <- heartsvg.out$score_col
    }
  }

  result
}

make_key <- function(pattern, peak.t.value, sd.r.value, replicate) {
  paste(pattern, as.numeric(peak.t.value), as.numeric(sd.r.value), as.integer(replicate), sep = "\r")
}

param.keys <- make_key(params$pattern, params$peak.t.value, params$sd.r.value, params$replicate)
outside.existing <- data.frame()
if(nrow(existing.results) > 0) {
  existing.keys <- make_key(
    existing.results$pattern,
    existing.results$peak.t.value,
    existing.results$sd.r.value,
    existing.results$replicate
  )
  outside.existing <- existing.results[!(existing.keys %in% param.keys), , drop = FALSE]
}

results <- vector("list", nrow(params))
for(i in seq_len(nrow(params))) {
  in.file <- input_file_for(
    pattern = params$pattern[i],
    peak.t.value = params$peak.t.value[i],
    sd.r.value = params$sd.r.value[i],
    replicate = params$replicate[i]
  )
  results[[i]] <- existing_result(
    pattern = params$pattern[i],
    peak.t.value = params$peak.t.value[i],
    sd.r.value = params$sd.r.value[i],
    replicate = params$replicate[i],
    in.file = in.file
  )
}

for(i in seq_len(nrow(params))) {
  results[[i]] <- run_one_dataset(
    pattern = params$pattern[i],
    peak.t.value = params$peak.t.value[i],
    sd.r.value = params$sd.r.value[i],
    replicate = params$replicate[i]
  )

  results.df <- do.call(rbind, results[seq_len(i)])
  if(i < nrow(params)) {
    results.df <- rbind(results.df, do.call(rbind, results[(i + 1):nrow(params)]))
  }
  if(nrow(outside.existing) > 0) {
    results.df <- rbind(outside.existing, results.df)
  }
  write.csv(results.df, results.file, row.names = FALSE)
  message("Checkpointed ", i, " / ", nrow(params), " rows to ", results.file)
}

message("Wrote final R model results to ", normalizePath(results.file))
