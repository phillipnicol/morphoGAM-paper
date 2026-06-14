
### Fit the the curve 

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

#t.smooth <- mgcv::gam(fit$xyt$t ~ s(fit$xyt$x, fit$xyt$y, bs="tp", k=100))

#fit$xyt$t <- predict(t.smooth)

p2 <- data.frame(x=fit$xyt$x,y=fit$xyt$y,color=fit$xyt$t) |>
  ggplot(aes(x=x,y=y,color=color)) + geom_point() +
  scale_color_gradientn(values=c(0,0.5,1),
                        colors=c("blue","grey90", "red"))+
  theme_bw() +
  ggtitle("Coordinate along curve") + labs(color="t")

#r <- ifelse(fit$xyt$r < 0.5, fit$xyt$r - 0.5, -1)
r <- fit$xyt$r 

p2 <- data.frame(x=fit$xyt$x,y=fit$xyt$y,color=r) |>
  ggplot(aes(x=x,y=y,color=color)) + geom_point() +
  scale_color_gradientn(values=c(0,0.5,1),
                        colors=c("blue","grey90", "red"))+
  theme_bw() +
  ggtitle("Distance from curve") + labs(color="Distance")





library(ggplot2)

df <- data.frame(x=fit$xyt$x, y=fit$xyt$y, r=fit$xyt$r)

p <- ggplot() +

  # Tumor (or graded) cells: r < 0.5
  geom_point(
    data = df[df$r > 0.5, ],
    aes(x = x, y = y, color = r),
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

library(ggpubr)
ggarrange(p1, p,nrow=1)

p1 <- fit$curve.plot + ggtitle("Fitted curve along tumor boundary")
p <- p + ggtitle("Second coordinate (distance from tumor edge)")

fit.full <- fit
save(fit.full, file="../data/sample_A_fit_full2.RData")



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

kxs <- which(fit$xyt$r > 0.5)

fit$xyt <- fit$xyt[kxs,]
mat_r <- mat_r[,kxs]


gene <- mat_r["hsp70l",]
data.frame(x=fit$xyt$x, y=fit$xyt$y, color=gene) |>
  ggplot(aes(x=x,y=y,color=color)) + geom_point() +
  scale_color_gradient(low="lightgrey", high="red") + theme_void()



rs <- rowSums(mat_r) 
mat_r <- mat_r[which(rs >= 25),]
rownames(mat_r) <- make.unique(rownames(mat_r))
mat_r <- as.matrix(mat_r)

#mgam <- MorphoGAM::MorphoGAM(Y=mat_r,
#                  curve.fit = fit,
#                  design = y ~ s(t, bs="cc") + s(r))


by.t <- ifelse(fit$xyt$r < 0.7, 1, 0)

data.frame(x=fit$xyt$x, y=fit$xyt$y, color=by.t) |>
  ggplot(aes(x=x,y=y,color=color)) + geom_point()  + theme_void()

mgam <- morphoGAM_test(Y=mat_r,
                  by.t=by.t,
                  curve.fit = fit,
                  design = y ~ s(t, bs="cc", by=by.t) + s(r))

save(mgam, file="../data/sample_A_mgam_results2.RData")

library(tidyverse)
library(RColorBrewer)
library(ggrepel)
#subset
qxs <- which(fit$xyt$r < 0.7)
mgam.sub <- mgam
mgam.sub$fxs.t <- mgam$fxs.t[,qxs]
mgam.sub$fxs.r <- mgam$fxs.r[,qxs]
mgam.sub$fpca.t <- irlba::irlba(mgam.sub$fxs.t,nv=5)
mgam.sub$fpca.r <- irlba::irlba(mgam.sub$fxs.r,nv=5)
fit.sub <- fit
fit.sub$xyt <- fit$xyt[qxs,]

p.fpc <- plotFPCloading_temp(mgam_object=mgam.sub, curve.fit=fit.sub, type="t", L=3, num_genes=10, sign="-")

library(gprofiler2)

# -----------------------------
# 1. Zebrafish invasive-edge genes
# -----------------------------

edge_genes <- rownames(mat_r)[order(mgam.sub$fpca.t$u[,3], decreasing = TRUE)[1:50]]
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


write.csv(
  as.matrix(human_results),
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

p.geneset <- ggplot(
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
      x = "Fraction of genes in term",
      y = NULL,
      size = "Overlap",
      color = "-log10(p-val)"
    )



p.geneset

library(ggpubr)

p.all <- ggarrange(ggarrange(p1, p, nrow=1), ggarrange(p.fpc, p.geneset,nrow=1), nrow=2, ncol=1)

mgam$results |> filter(peak_r_location > 0.7) |> arrange(desc(peak.r)) |> head(n=20)

gene <- "hsp70l"
plot(fit$xyt$r, mat_r[gene,])
plot(fit$xyt$r, mgam$fxs.r[gene,])
plot(fit$xyt$t, mat_r[gene,])
plot(fit$xyt$t, mgam$fxs.t[gene,])
p.hsp.raw <- data.frame(x=fit$xyt$x, y=fit$xyt$y, color=mat_r[gene,]) |>
  ggplot(aes(x=x,y=y,color=color)) + geom_point() +
  scale_color_gradient(low="lightgrey", high="red") + theme_void() + 
  ggtitle(paste0(gene, " expression")) + labs(color="Raw count")

p.hsp.r <- MorphoGAM:::plotGAMestimates(mat_r,
genes=c("hsp70l"),
curve_fit=fit,
mgam=mgam,
type="r"
)

p.all <- ggarrange(ggarrange(p1, p, nrow=1, labels=c("a","b")), 
                    ggarrange(p.fpc, p.geneset, labels=c("c","d"), nrow=1),
                    ggarrange(p.hsp.raw, p.hsp.r, nrow=1, labels=c("e","f")),
                    nrow=3,ncol=1, heights=c(1.25, 1, 1)
)

ggsave(p.all,
       filename="../plots/sample_A_figure.png",
       width=12, height=12,
       units="in")
