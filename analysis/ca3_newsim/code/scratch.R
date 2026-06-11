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

#save xy as csv
write.csv(xy, file = "../data/ca3_xy.csv", row.names = FALSE)


G <- 1000 #Number of genes

library(splines)

df.t <- 20
T.basis <- ns(fit$xyt$t, df=df.t)
R.basis <- ns(fit$xyt$r, df=3)

n <- ncol(Y)
p <- 250

sd.t = 1
sd.r = 1
theta <- 10


Mu <- matrix(0, nrow = p, ncol = n)
Y.new <- Mu

#Sample 1 hot vector of length ncol(T.basis)
v <- rep(0, ncol(T.basis))
v[sample(2:(ncol(T.basis)-1), 1)] <- 1

alpha1 <- -1
peak.t <- 1.5
Mu[1,] <- alpha1 + T.basis %*% (peak.t*v)
Y.new[1,] <- rnbinom(n, mu = exp(Mu[1,]), size = theta)

for(j in 2:p) {
  print(j)
  Mu[j,] <- R.basis %*% rnorm(ncol(R.basis),sd=sd.r)

  #Generate Y.new[j,] from Nbin

  alphaj <- rnorm(1)
  slopej <- rnorm(1,sd=5)

  Mu[j,] <- alphaj + slopej*fit$xyt$r
  
  Y.new[j,] <- rnbinom(n, mu = exp(Mu[j,]), size = theta)
}

Y.new <- as.matrix(Y.new)
rownames(Y.new) <- 1:p

#Run the MorphoGAM

mgam <- MorphoGAM::MorphoGAM(Y=Y.new,
            curve.fit = fit,
            design = y ~ s(t,bs="cr",k=20) + s(r),
            return.fx = TRUE)

mgam$results |> arrange(desc(peak.t)) |> head(n=5)

#mgam$results |> head()

#Run the SPARK 
library(SPARK)
spark.res <- SPARK::sparkx(count_in=Y.new, 
                         locus_in = xy)

spark.res$res_mtest |> arrange(adjustedPval) |> head(n=5)

#Ranking of 1 
which(rownames(spark.res$res_mtest |> arrange(adjustedPval)) == "1")

#Run the nnSVG

library(nnSVG)
logCPM <- log(sweep(Y.new, MARGIN=2, STATS = colSums(Y.new), FUN="/")+1)
nn.svg <- nnSVG(input=logCPM, spatial_coords = xy,
                verbose=TRUE)

nn.svg |> arrange(pdj) |> head(n=5)


#Run Another SVG method 






### Testing datasets created in create_datasets.R
Y_sim <- read.csv("../data/ca3_simulated_dataset_peakt_1_sd.r_2_rep_01.csv")
Y_sim <- as.matrix(Y_sim[,-1])
rownames(Y_sim) <- 1:nrow(Y_sim)

mgam <- MorphoGAM::MorphoGAM(Y=Y_sim,
            curve.fit = fit,
            design = y ~ s(t,bs="cr",k=20) + s(r),
            return.fx = TRUE)

mgam$results |> arrange(desc(peak.t)) |> head(n=10)



#Run the SPARK 
library(SPARK)
spark.res <- SPARK::sparkx(count_in=Y_sim, 
                         locus_in = xy)

spark.res$res_mtest |> arrange(adjustedPval) |> head(n=5)

#Ranking of 1 
which(rownames(spark.res$res_mtest |> arrange(adjustedPval)) == "1")




#Run the nnSVG

library(nnSVG)
logCPM <- log(sweep(Y.new, MARGIN=2, STATS = colSums(Y_sim), FUN="/")+1)
nn.svg <- nnSVG(input=logCPM, spatial_coords = xy,
                verbose=TRUE)

nn.svg[1,12]

#Save t coord
write.csv(fit$xyt$t, file = "../data/ca3_t_coords.csv", row.names = FALSE)
write.csv(fit$xyt$r, file = "../data/ca3_r_coords.csv", row.names = FALSE)


### Help spateo get points
which.min(fit$xyt$r[which(fit$xyt$t < 0.01)])
which.max(fit$xyt$r[which(fit$xyt$t < 0.01)])
which.min(fit$xyt$r[which(fit$xyt$t > 0.99)])
which.max(fit$xyt$r[which(fit$xyt$t > 0.99)])

fit$xyt$spot <- 1:nrow(fit$xyt)
fit$xyt |> filter(t < 0.05) |> arrange(desc(r))
fit$xyt |> filter(t < 0.05) |> arrange(r)
fit$xyt |> filter(t > 0.95) |> arrange(desc(r))
fit$xyt |> filter(t > 0.95) |> arrange(r)




### Testing datasets created in create_datasets.R
Y_sim <- read.csv("../data/ca3_simulation_gradient/ca3_simulated_dataset_peakt_1_sd.r_2_rep_02.csv")
Y_sim <- as.matrix(Y_sim[,-1])
rownames(Y_sim) <- 1:nrow(Y_sim)

#Run the SPARK 
library(SPARK)
spark.res <- SPARK::sparkx(count_in=Y_sim, 
                         locus_in = xy)

spark.res$res_mtest |> arrange(adjustedPval) |> head(n=5)
which(rownames(spark.res$res_mtest |> arrange(adjustedPval)) == "1")

mgam <- MorphoGAM::MorphoGAM(Y=Y_sim,
            curve.fit = fit,
            design = y ~ s(t,bs="cr",k=20) + s(r),
            return.fx = TRUE)

mgam$results |> arrange(desc(peak.t)) |> head(n=10)

linear.functional <- apply(mgam$fxs.t, 1, function(x) x[which.min(fit$xyt$t)] - x[which.max(fit$xyt$t)])
