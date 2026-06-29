#load("../data/granule_curve.RData")
#out <- fit

setwd(here::here("analysis/granule_t_sim/code"))


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

  #Load ground truth
  load(here::here("analysis", "principal_curves", "data", "granule_ground_truth.Rda"))
  fit <- ground_truth
  #fit <- CurveFinder(xy, knn = 9)

  dir.create(file.path("..", "data"), showWarnings = FALSE, recursive = TRUE)
  write.csv(as.data.frame(xy), file = "../data/granule_xy.csv", row.names = FALSE)
  write.csv(fit$xyt$t, file = "../data/granule_t_coords.csv", row.names = FALSE)
  write.csv(fit$xyt$r, file = "../data/granule_r_coords.csv", row.names = FALSE)

  digitized <- data.frame(
    xcoord = xy[, 1],
    ycoord = xy[, 2],
    spateo_longitudinal = fit$xyt$t,
    spateo_radial = fit$xyt$r,
    spateo_longitudinal_raw = fit$xyt$t,
    spateo_radial_raw = fit$xyt$r
  )
  write.csv(digitized, file = "../data/granule_xy_digitized.csv", row.names = TRUE)

  list(xy = xy, fit = fit)
}

granule <- prepare_granule_curve()
xy <- granule$xy
fit <- granule$fit
out <- fit


set.seed(1)

### Part I: Create spike gene datasets 

M <- 10 #Number of simulated datasets
peak.t <- c(1, 1.5, 2)
p <- 1000 #Number of simulated genes per dataset
sd.r <- c(0.5, 1, 2)


library(splines)

df.t <- 10
T.basis <- ns(fit$xyt$t, df=df.t)

n <- nrow(xy)

theta <- 10


out.dir <- file.path("..", "data", "granule_simulation_peak")
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)

simulate_counts <- function(peak.t.value, sd.r.value) {
  Mu <- matrix(0, nrow = p, ncol = n)
  Y.new <- Mu

  # Sample 1 hot vector of length ncol(T.basis)
  v <- rep(0, ncol(T.basis))
  v[sample(2:(ncol(T.basis) - 1), 1)] <- 1

  alpha1 <- rnorm(1, mean=-1)
  Mu[1,] <- alpha1 + T.basis %*% (peak.t.value * v)
  Y.new[1,] <- rnbinom(n, mu = exp(Mu[1,]), size = theta)

  for(j in 2:p) {
    # Generate Y.new[j,] from Nbin
    alphaj <- rnorm(1)
    slopej <- rnorm(1, sd = sd.r.value) #Change this as needed

    Mu[j,] <- alphaj + slopej * fit$xyt$r
    Y.new[j,] <- rnbinom(n, mu = exp(Mu[j,]), size = theta)
  }

  Y.new <- as.matrix(Y.new)
  rownames(Y.new) <- 1:p
  Y.new
}

output.files <- character()

params <- expand.grid(peak.t.value = peak.t, sd.r.value = sd.r)

for(i in 1:nrow(params)) {
  peak.t.value <- params[i, 1]
  sd.r.value <- params[i, 2]
  peak.t.label <- gsub("\\.", "p", format(peak.t.value, trim = TRUE, scientific = FALSE))
  sd.r.label <- gsub("\\.", "p", format(sd.r.value, trim = TRUE, scientific = FALSE))

  for(m in seq_len(M)) {
    message("Simulating peak.t = ", peak.t.value, ", replicate = ", m)

    Y.new <- simulate_counts(peak.t.value, sd.r.value)
    out.file <- file.path(
      out.dir,
      sprintf("granule_simulated_dataset_peakt_%s_sd.r_%s_rep_%02d.csv", peak.t.label, sd.r.label, m)
    )

    write.csv(Y.new, out.file)
    output.files <- c(output.files, out.file)
  }
}

message("Wrote ", length(output.files), " datasets to ", normalizePath(out.dir))



## Part II : Linear slope in t 


set.seed(1)


M <- 10 #Number of simulated datasets
peak.t <- c(1, 1.5, 2)
p <- 1000 #Number of simulated genes per dataset
sd.r <- c(0.5, 1, 2)


library(splines)

df.t <- 10
T.basis <- ns(fit$xyt$t, df=df.t)

n <- nrow(xy)

theta <- 10


out.dir <- file.path("..", "data", "granule_simulation_gradient")
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)

simulate_counts <- function(peak.t.value, sd.r.value) {
  Mu <- matrix(0, nrow = p, ncol = n)
  Y.new <- Mu

  # Sample 1 hot vector of length ncol(T.basis)
  v <- rep(0, ncol(T.basis))
  v[sample(2:(ncol(T.basis) - 1), 1)] <- 1

  alpha1 <- rnorm(1, mean=-1)
  #Mu[1,] <- alpha1 + T.basis %*% (peak.t.value * v)
  Mu[1,] <- alpha1 + peak.t.value*fit$xyt$t
  Y.new[1,] <- rnbinom(n, mu = exp(Mu[1,]), size = theta)

  for(j in 2:p) {
    # Generate Y.new[j,] from Nbin
    alphaj <- rnorm(1, mean=-1)
    slopej <- rnorm(1, sd = sd.r.value) #Change this as needed

    Mu[j,] <- alphaj + slopej * fit$xyt$r
    Y.new[j,] <- rnbinom(n, mu = exp(Mu[j,]), size = theta)
  }

  Y.new <- as.matrix(Y.new)
  rownames(Y.new) <- 1:p
  Y.new
}

output.files <- character()

params <- expand.grid(peak.t.value = peak.t, sd.r.value = sd.r)

for(i in 1:nrow(params)) {
  peak.t.value <- params[i, 1]
  sd.r.value <- params[i, 2]
  peak.t.label <- gsub("\\.", "p", format(peak.t.value, trim = TRUE, scientific = FALSE))
  sd.r.label <- gsub("\\.", "p", format(sd.r.value, trim = TRUE, scientific = FALSE))

  for(m in seq_len(M)) {
    message("Simulating peak.t = ", peak.t.value, ", replicate = ", m)

    Y.new <- simulate_counts(peak.t.value, sd.r.value)
    out.file <- file.path(
      out.dir,
      sprintf("granule_simulated_dataset_peakt_%s_sd.r_%s_rep_%02d.csv", peak.t.label, sd.r.label, m)
    )

    write.csv(Y.new, out.file)
    output.files <- c(output.files, out.file)
  }
}

message("Wrote ", length(output.files), " datasets to ", normalizePath(out.dir))



## Part III: Non-parametric


set.seed(1)


M <- 10 #Number of simulated datasets
peak.t <- c(1, 1.5, 2)
p <- 1000 #Number of simulated genes per dataset
sd.r <- c(0.5, 1, 2)


library(splines)

df.t <- 10
T.basis <- ns(fit$xyt$t, df=df.t)

n <- nrow(xy)

theta <- 10


out.dir <- file.path("..", "data", "granule_simulation_nonparam")
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)

simulate_counts <- function(peak.t.value, sd.r.value) {
  Mu <- matrix(0, nrow = p, ncol = n)
  Y.new <- Mu

  # Sample 1 hot vector of length ncol(T.basis)
  v <- rep(0, ncol(T.basis))
  v[sample(2:(ncol(T.basis) - 1), 1)] <- 1

  alpha1 <- rnorm(1, mean=-1)
  #Mu[1,] <- alpha1 + T.basis %*% (peak.t.value * v)
  Mu[1,] <- alpha1 + T.basis %*% rnorm(ncol(T.basis))
  Mu[1,] <- Mu[1,]/max(abs(Mu[1,]))*peak.t.value #Scale to have desired peak.t.value

  Y.new[1,] <- rnbinom(n, mu = exp(Mu[1,]), size = theta)

  for(j in 2:p) {
    # Generate Y.new[j,] from Nbin
    alphaj <- rnorm(1, mean=-1)
    slopej <- rnorm(1, sd = sd.r.value) #Change this as needed

    Mu[j,] <- alphaj + slopej * fit$xyt$r
    Y.new[j,] <- rnbinom(n, mu = exp(Mu[j,]), size = theta)
  }

  Y.new <- as.matrix(Y.new)
  rownames(Y.new) <- 1:p
  Y.new
}

output.files <- character()

params <- expand.grid(peak.t.value = peak.t, sd.r.value = sd.r)

for(i in 1:nrow(params)) {
  peak.t.value <- params[i, 1]
  sd.r.value <- params[i, 2]
  peak.t.label <- gsub("\\.", "p", format(peak.t.value, trim = TRUE, scientific = FALSE))
  sd.r.label <- gsub("\\.", "p", format(sd.r.value, trim = TRUE, scientific = FALSE))

  for(m in seq_len(M)) {
    message("Simulating peak.t = ", peak.t.value, ", replicate = ", m)

    Y.new <- simulate_counts(peak.t.value, sd.r.value)
    out.file <- file.path(
      out.dir,
      sprintf("granule_simulated_dataset_peakt_%s_sd.r_%s_rep_%02d.csv", peak.t.label, sd.r.label, m)
    )

    write.csv(Y.new, out.file)
    output.files <- c(output.files, out.file)
  }
}

message("Wrote ", length(output.files), " datasets to ", normalizePath(out.dir))
