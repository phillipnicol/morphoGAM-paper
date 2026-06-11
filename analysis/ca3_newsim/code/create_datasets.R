#load("../data/ca3_curve.RData")
#out <- fit

setwd(here::here("analysis/ca3_newsim/code"))


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

### Part I: Create spike gene datasets 

M <- 10 #Number of simulated datasets
peak.t <- c(1, 1.5, 2)
p <- 1000 #Number of simulated genes per dataset
sd.r <- c(0.5, 1, 2)


library(splines)

df.t <- 10
T.basis <- ns(fit$xyt$t, df=df.t)

n <- ncol(Y)

theta <- 10


out.dir <- file.path("..", "data", "ca3_simulation_peak")
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)

simulate_counts <- function(peak.t.value, sd.r.value) {
  Mu <- matrix(0, nrow = p, ncol = n)
  Y.new <- Mu

  # Sample 1 hot vector of length ncol(T.basis)
  v <- rep(0, ncol(T.basis))
  v[sample(2:(ncol(T.basis) - 1), 1)] <- 1

  alpha1 <- -1
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
      sprintf("ca3_simulated_dataset_peakt_%s_sd.r_%s_rep_%02d.csv", peak.t.label, sd.r.label, m)
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

n <- ncol(Y)

theta <- 10


out.dir <- file.path("..", "data", "ca3_simulation_gradient")
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
      sprintf("ca3_simulated_dataset_peakt_%s_sd.r_%s_rep_%02d.csv", peak.t.label, sd.r.label, m)
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

n <- ncol(Y)

theta <- 10


out.dir <- file.path("..", "data", "ca3_simulation_nonparam")
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
      sprintf("ca3_simulated_dataset_peakt_%s_sd.r_%s_rep_%02d.csv", peak.t.label, sd.r.label, m)
    )

    write.csv(Y.new, out.file)
    output.files <- c(output.files, out.file)
  }
}

message("Wrote ", length(output.files), " datasets to ", normalizePath(out.dir))
