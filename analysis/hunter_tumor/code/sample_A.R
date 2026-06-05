
setwd(here::here("analysis/hunter_tumor/code"))

#source("../../../R/curve_finder_i_try.R")
#source("../../../R/helpers.R")

library(ggplot2)

xy_raw <- read.csv("../data/coords_sampleA.csv", header = FALSE)

full_w <- 31870
full_h <- 32588

small_w <- 785
small_h <- 768

image_path <- "../data/GSM4838131_Visium_A_image.tif"

ixs <- which(xy_raw$V2 == 1)

x_full <- xy_raw[ixs, 5]
y_full <- xy_raw[ixs, 6]

# aspect-preserving shrink matched to height
s <- small_h / full_h

expected_w <- full_w * s
pad_x <- (small_w - expected_w) / 2

# manual correction: dots looked too far right
x_shift <- -10
y_shift <- 7

xy_small <- cbind(
  x = x_full * s + pad_x + x_shift,
  y = y_full * s + y_shift
)

# Do NOT flip y here, because CurveFinderInteractive2 does it
fit <- MorphoGAM:::CurveFinderInteractive2(
  xy = xy_small,
  loop = TRUE,
  image_path = image_path,
  flip_y = FALSE
)





p1 <- fit$curve.plot + ggtitle("Fitted curve along tumor boundary")

t.smooth <- mgcv::gam(fit$xyt$t ~ s(fit$xyt$x, fit$xyt$y, bs="tp", k=100))

fit$xyt$t <- predict(t.smooth)

p2 <- data.frame(x=fit$xyt$x,y=fit$xyt$y,color=fit$xyt$t) |>
  ggplot(aes(x=x,y=y,color=color)) + geom_point() +
  scale_color_gradientn(values=c(0,0.5,1),
                        colors=c("blue","grey90", "red"))+
  theme_bw() +
  ggtitle("Coordinate along curve") + labs(color="t")

r <- ifelse(fit$xyt$r < 0.5, fit$xyt$r - 0.5, -1)

p2 <- data.frame(x=fit$xyt$x,y=fit$xyt$y,color=r) |>
  ggplot(aes(x=x,y=y,color=color)) + geom_point() +
  scale_color_gradientn(values=c(0,0.5,1),
                        colors=c("blue","grey90", "red"))+
  theme_bw() +
  ggtitle("Distance from curve") + labs(color="Distance")





library(ggplot2)

df <- data.frame(x=fit$xyt$x, y=-fit$xyt$y, r=fit$xyt$r)

p <- ggplot() +

  # Tumor (or graded) cells: r < 0.5
  geom_point(
    data = df[df$r > 0.5, ],
    aes(x = x, y = y, color = 0.5 - r),
    size = 1
  ) +

  # Non-tumor cells: r > 0.5
  geom_point(
    data = df[df$r < 0.5, ],
    aes(x = x, y = y),
    color = "grey90",
    size = 1
  ) +

  scale_color_gradient(
    low = "darkblue",
    high = "lightblue",
    name = "r"
  ) +

  guides(
    color = guide_colorbar(order = 1),
    fill  = "none"
  ) +

  labs(
    color = "r",
    caption = "Grey points: non-tumor cell"
  ) +

  theme_bw()



p1 <- fit$curve.plot + ggtitle("Fitted curve along tumor boundary")
p <- p + ggtitle("Second coordinate (distance from tumor edge)")

kxs <- which(fit$xyt$r > 0.5)


library(anndata)

library(rhdf5)

library(rhdf5)
library(Matrix)

file <- "../data/GSM4838132_Visium_Sample_B_filtered_feature_bc_matrix.h5"

data    <- h5read(file, "matrix/data")
indices <- h5read(file, "matrix/indices")
indptr  <- h5read(file, "matrix/indptr")
shape   <- h5read(file, "matrix/shape")

# Build sparse matrix (note: +1 because R is 1-based)
mat <- sparseMatrix(
  i = indices + 1,
  p = indptr,
  x = data,
  dims = shape
)

genes <- h5read(file, "matrix/features/name")
cells <- h5read(file, "matrix/barcodes")

rownames(mat) <- genes
colnames(mat) <- cells


cell_meta <- data.frame(
  barcode = h5read(file, "matrix/barcodes"),
  stringsAsFactors = FALSE
)

head(cell_meta)


mat_r <- as(mat, "RsparseMatrix")

jxs <- match(xy_raw[ixs,]$V1, colnames(mat_r))

mat_r <- mat_r[,jxs]

gene <- mat_r[which(rownames(mat_r) == "RPL41"),]

fit$xyt <- fit$xyt[kxs,]
mat_r <- mat_r[,kxs]


data.frame(x=xy_small[,1], y=xy_small[,2], color=gene) |>
  ggplot(aes(x=x,y=y,color=color)) + geom_point() +
  scale_color_gradient(low="lightgrey", high="red") + theme_void()



rs <- rowSums(mat_r) 
mat_r <- mat_r[which(rs >= 25),]
rownames(mat_r) <- make.unique(rownames(mat_r))
mat_r <- as.matrix(mat_r)

mgam <- MorphoGAM::MorphoGAM(Y=mat_r,
                  curve.fit = fit,
                  design = y ~ s(t, bs="cc") + s(r))


by.t <- ifelse(fit$xyt$r < 0.7, 1, 0)

data.frame(x=fit$xyt$x, y=fit$xyt$y, color=by.t) |>
  ggplot(aes(x=x,y=y,color=color)) + geom_point()  + theme_void()

mgam <- morphoGAM_test(Y=mat_r,
                  by.t=by.t,
                  curve.fit = fit,
                  design = y ~ s(t, bs="cc", by=by.t) + s(r))


plot(fit$xyt$r, mat_r[4,])

peak_r_location <- fit$xyt$r[apply(mgam$fxs.r, 1, which.max)]

rownames(mat_r)[order(abs(peak_r_location - 0.5), decreasing = FALSE)[1:10]]
names(peak_r_location) <- rownames(mat_r)
which(peak_r_location < 0.6 & mgam$results$peak.r > log(2))

plot(fit$xyt$r, mat_r["icn",])

plot(fit$xyt$t, mat_r["vmhc",])

peak_t_location <- fit$xyt$t[apply(mgam$fxs.t, 1, which.max)]
names(peak_t_location) <- rownames(mat_r)

which(peak_t_location > 0.2 & peak_t_location < 0.7 & mgam$results$peak.t > log(2) &
  peak_r_location < 0.6 & mgam$results$peak.r > log(2))

which(peak_r_location > 0.75 & mgam$results$peak.r > log(2))

data.frame(x=fit$xyt$x, y=fit$xyt$y, color=mat_r["nkx3-1",]) |>
  ggplot(aes(x=x,y=y,color=color)) + geom_point() +
  scale_color_gradient(low="lightgrey", high="red") + theme_void()

plot(fit$xyt$t, mgam$fxs.t["rbp4",])
plot(fit$xyt$r, mgam$fxs.r["rbp4",])

plot(fit$xyt$t, mat_r["pvalb4",])

which(peak_t_location > 0.2 & peak_t_location < 0.7 & mgam$results$peak.t > log(2))

MorphoGAM:::plotFPCloading(mgam, curve.fit=fit, type="t", L=2, num_genes=4)

MorphoGAM:::plotFPCloading(mgam, curve.fit=fit, type="r", L=1, num_genes=10)

plotFPCloading_temp(mgam, curve.fit=fit, type="t", L=3, num_genes=10, sign="-")


plotFPCloading_temp(mgam, curve.fit=fit, type="t", L=1, num_genes=10, sign="+")
plotFPCloading_temp(mgam, curve.fit=fit, type="t", L=5, num_genes=10, sign="-")


plotFPCloading_temp(mgam, curve.fit=fit, type="r", L=3, num_genes=10, sign="+")

data.frame(x=fit$xyt$x, y=fit$xyt$y, color=mat_r["col1a1b",]) |>
  ggplot(aes(x=x,y=y,color=color)) + geom_point() +
  scale_color_gradient(low="lightgrey", high="red") + theme_void()

gene <- "ano11"
plot(fit$xyt$r, mat_r[gene,])
plot(fit$xyt$r, mgam$fxs.r[gene,])
plot(fit$xyt$t, mat_r[gene,])
plot(fit$xyt$t, mgam$fxs.t[gene,])
data.frame(x=fit$xyt$x, y=fit$xyt$y, color=mat_r[gene,]) |>
  ggplot(aes(x=x,y=y,color=color)) + geom_point() +
  scale_color_gradient(low="lightgrey", high="red") + theme_void()

max.ix <- which.max(fit$xyt$r)
min.ix <- which.min(fit$xyt$r)
diff.stat <- apply(mgam$fxs.r, 1, function(fx) {
  fx[max.ix] - fx[min.ix]
})

diff.stat.exp <- exp(mgam$results$intercept + diff.stat)



fpca.old <- mgam$fpca.t



mgam$fpca.t <- irlba::irlba(mgam$fxs.t[,by.t==1],nv=5)



mgam$results |> arrange(desc(peak.t)) |> head(n=20)


which(peak_r_location > 0.7 & mgam$results$peak.r > 0.5*log(2))

gene <- "spp1"
plot(fit$xyt$r, mat_r[gene,])
plot(fit$xyt$r, mgam$fxs.r[gene,])
plot(fit$xyt$t, mat_r[gene,])
plot(fit$xyt$t, mgam$fxs.t[gene,])
data.frame(x=fit$xyt$x, y=fit$xyt$y, color=mat_r[gene,]) |>
  ggplot(aes(x=x,y=y,color=color)) + geom_point() +
  scale_color_gradient(low="lightgrey", high="red") + theme_void()



install.packages("gprofiler2")
library(gprofiler2)

genes <- which(peak_r_location > 0.7 & mgam$results$peak.r > 0.3*log(2)) |> names()


genes <- mgam$results |> filter(peak_r_location > 0.7) |>
  arrange(desc(peak.r)) |> head(n=50) |> rownames()

# Direct zebrafish enrichment
res <- gost(
  query = genes,
  organism = "drerio",
  sources = c("GO:BP", "GO:CC", "REAC", "WP"),
  significant = FALSE
)

res$result[1:10, c("source", "term_name", "p_value", "intersection")]

gconvert(
  query = genes,
  organism = "drerio"
)









#### Subset


genes <- rownames(mat_r)[order(mgam.sub$fpca.t$u[,3], decreasing = FALSE)[1:50]]
# Direct zebrafish enrichment
res <- gost(
  query = genes,
  organism = "drerio",
  sources = c("GO:BP", "GO:CC", "REAC", "WP"),
  significant = FALSE
)

res$result[1:10, c("source", "term_name", "p_value", "intersection")]

qxs <- which(fit$xyt$r < 0.7)
gene <- "col5a1"
plot(fit$xyt$t[qxs], mgam$fxs.t[gene,qxs])

#subset
mgam.sub <- mgam
mgam.sub$fxs.t <- mgam$fxs.t[,qxs]
mgam.sub$fxs.r <- mgam$fxs.r[,qxs]
mgam.sub$fpca.t <- irlba::irlba(mgam.sub$fxs.t,nv=5)
mgam.sub$fpca.r <- irlba::irlba(mgam.sub$fxs.r,nv=5)
fit.sub <- fit
fit.sub$xyt <- fit$xyt[qxs,]

plotFPCloading_temp(mgam_object=mgam.sub, curve.fit=fit.sub, type="t", L=3, num_genes=10, sign="-")







# Install if needed
install.packages("gprofiler2")

library(gprofiler2)

# -----------------------------
# 1. Zebrafish invasive-edge genes
# -----------------------------

edge_genes <- genes
# Remove duplicated symbols
edge_genes <- unique(edge_genes)

# -----------------------------
# 2. Check zebrafish gene mapping
# -----------------------------

mapped <- gconvert(
  query = edge_genes,
  organism = "drerio"
)

mapped

# See which genes failed to map
unmapped <- setdiff(edge_genes, mapped$input)
unmapped

# Use only successfully mapped genes
edge_mapped <- unique(mapped$input)

# -----------------------------
# 3. Direct zebrafish enrichment
# -----------------------------

zf_enrich <- gost(
  query = edge_mapped,
  organism = "drerio",
  sources = c("GO:BP", "GO:CC", "GO:MF", "REAC", "WP"),
  significant = FALSE,
  user_threshold = 1
)

zf_results <- zf_enrich$result

# View top results
zf_results[1:30, c(
  "source",
  "term_name",
  "p_value",
  "term_size",
  "query_size",
  "intersection"
)]

# Save results
write.csv(
  zf_results,
  "../data/zebrafish_edge_enrichment_results.csv",
  row.names = FALSE
)

# -----------------------------
# 4. Human ortholog conversion
# -----------------------------

orth <- gorth(
  query = edge_mapped,
  source_organism = "drerio",
  target_organism = "hsapiens"
)

orth[, c("input", "ortholog_name")]

human_genes <- unique(orth$ortholog_name)
human_genes <- human_genes[!is.na(human_genes)]

human_genes

# -----------------------------
# 5. Human ortholog enrichment
# -----------------------------

human_enrich <- gost(
  query = human_genes,
  organism = "hsapiens",
  sources = c("GO:BP", "GO:CC", "GO:MF", "REAC", "WP", "HPA"),
  significant = FALSE,
  user_threshold = 1
)

human_results <- human_enrich$result

human_results <- as.matrix(human_results)

write.csv(
  human_results,
  "../data/human_ortholog_edge_enrichment_results.csv",
  row.names = FALSE
)

# -----------------------------
# 6. Filter for biologically relevant terms
# -----------------------------

keywords <- c(
  "extracellular",
  "matrix",
  "collagen",
  "migration",
  "motility",
  "wound",
  "inflammatory",
  "immune",
  "macrophage",
  "metalloproteinase",
  "adhesion",
  "fibronectin",
  "basement",
  "invasion",
  "tissue remodeling"
)

pattern <- paste(keywords, collapse = "|")

relevant_zf <- zf_results[
  grepl(pattern, zf_results$term_name, ignore.case = TRUE),
]

relevant_human <- human_results[
  grepl(pattern, human_results$term_name, ignore.case = TRUE),
]

relevant_zf[, c("source", "term_name", "p_value", "intersection")]

relevant_human[, c("source", "term_name", "p_value", "intersection")]

# -----------------------------
# 7. Plot top human ortholog results
# -----------------------------

top <- human_results[order(human_results$p_value), ][1:20, ]

barplot(
  -log10(top$p_value),
  names.arg = top$term_name,
  las = 2,
  cex.names = 0.6,
  ylab = "-log10(p-value)",
  main = "Invasive-edge gene enrichment"
)



library(ggplot2)
library(dplyr)
plot_terms <- human_results %>%
  filter(
    grepl(
      "matrix|ECM|collagen|basement|extracellular|proteoglycan",
      term_name,
      ignore.case = TRUE
    )
  ) %>%
  arrange(p_value) %>%
  slice_head(n = 10) %>%
  mutate(
    term_name = factor(term_name,
                       levels = rev(term_name))
  )

ggplot(
  plot_terms,
  aes(
    x = intersection_size / query_size,
    y = term_name
  )
) +
  geom_point(
    aes(
      size = intersection_size,
      color = -log10(p_value)
    )
  ) +
  theme_classic()





library(dplyr)
library(ggplot2)

plot_terms <- human_results %>%
  filter(
    grepl(
      "matrix|ECM|collagen|basement|extracellular|proteoglycan",
      term_name,
      ignore.case = TRUE
    )
  ) %>%
  mutate(
    gene_ratio = intersection_size / query_size
  ) %>%
  arrange(desc(gene_ratio)) %>%
  slice_head(n = 15) %>%
  mutate(
    term_name = factor(
      term_name,
      levels = rev(term_name)
    )
  )

ggplot(
  plot_terms,
  aes(
    x = gene_ratio,
    y = term_name
  )
) +
  geom_point(
    aes(
      size = intersection_size,
      color = -log10(p_value)
    )
  ) +
  theme_classic() +
  labs(
    x = "Fraction of invasive-edge genes",
    y = NULL,
    size = "Overlap",
    color = "-log10(P)"
  )