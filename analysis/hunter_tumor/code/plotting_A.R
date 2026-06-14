# Load fit and mgam
setwd(here::here("analysis/hunter_tumor/code"))

load("../data/sample_A_fit_full2.RData")
load("../data/sample_A_mgam_results2.RData")



library(ggplot2)
library(grid)

panel_text_theme <- theme(
  plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
  axis.title = element_text(size = 12),
  axis.text = element_text(size = 11),
  legend.title = element_text(size = 11),
  legend.text = element_text(size = 10),
  plot.margin = margin(8, 8, 8, 8)
)

spatial_text_theme <- panel_text_theme +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.caption = element_text(size = 10, hjust = 0)
  )

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


#Curve plot

library(ggplot2)
fit <- fit.full

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


p1 <- fit$curve.plot +
  ggtitle("Fitted curve along tumor boundary") +
  guides(color = guide_colorbar(barheight = unit(35, "pt"), barwidth = unit(5, "pt"))) +
  spatial_text_theme
p <- p +
  ggtitle("r (distance from tumor edge)") +
  guides(color = guide_colorbar(barheight = unit(35, "pt"), barwidth = unit(5, "pt"))) +
  spatial_text_theme

library(ggpubr)
p.tumor <- ggarrange(p1, p,nrow=1)



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
if (anyNA(jxs)) {
  stop(
    "Some coordinate barcodes were not found in the count matrix: ",
    paste(head(xy_raw[ixs,]$V1[is.na(jxs)], 10), collapse = ", ")
  )
}

mat_r <- mat_r[,jxs]

kxs <- which(fit$xyt$r > 0.5)

fit$xyt <- fit$xyt[kxs,]
mat_r <- mat_r[,kxs]


rs <- rowSums(mat_r) 
mat_r <- mat_r[which(rs >= 25),]
rownames(mat_r) <- make.unique(rownames(mat_r))
mat_r <- as.matrix(mat_r)



mgam_gene_order <- rownames(mgam$fxs.r)
if (!identical(rownames(mat_r), mgam_gene_order)) {
  missing_genes <- setdiff(mgam_gene_order, rownames(mat_r))
  extra_genes <- setdiff(rownames(mat_r), mgam_gene_order)

  if (length(missing_genes) > 0 || length(extra_genes) > 0) {
    stop(
      "mat_r and mgam gene names do not match. ",
      "Missing from mat_r: ", paste(head(missing_genes, 10), collapse = ", "),
      "; extra in mat_r: ", paste(head(extra_genes, 10), collapse = ", ")
    )
  }

  mat_r <- mat_r[mgam_gene_order, , drop = FALSE]
}

stopifnot(
  identical(rownames(mat_r), rownames(mgam$fxs.r)),
  identical(rownames(mat_r), rownames(mgam$fxs.t)),
  ncol(mat_r) == nrow(fit$xyt),
  ncol(mgam$fxs.r) == nrow(fit$xyt),
  ncol(mgam$fxs.t) == nrow(fit$xyt)
)

#mgam <- MorphoGAM::MorphoGAM(Y=mat_r,
#                  curve.fit = fit,
#                  design = y ~ s(t, bs="cc") + s(r))


set.seed(1)
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
p.fpc <- p.fpc +
  ggtitle("FPC loading genes") +
  coord_cartesian(clip = "off") +
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.22))) +
  panel_text_theme +
  theme(
    plot.margin = margin(8, 28, 8, 8),
    axis.title.y = element_text(margin = margin(r = 6))
  )

library(gprofiler2)

peak_r_location <- fit$xyt$r[apply(mgam$fxs.r, 1, which.max)]


# -----------------------------
# 1. Zebrafish invasive-edge genes
# -----------------------------

edge_genes <- rownames(mgam.sub$fxs.t)[order(mgam.sub$fpca.t$u[,3], decreasing = FALSE)[1:50]]
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
    slice_head(n = 10) %>% #Change back to 15 if need
    mutate(
      term_label = stringr::str_wrap(term_name, width = 32),
      term_name = factor(
        term_name,
        levels = rev(term_name)
      ),
      term_label = factor(
        term_label,
        levels = rev(unique(term_label))
      )
    )

p.geneset <- ggplot(
    plot_terms,
    aes(
      x = gene_ratio,
      y = term_label
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
    ) +
    guides(
      size = guide_legend(order = 1, title.position = "top"),
      color = guide_colorbar(order = 2, title.position = "top",
                             barheight = unit(38, "pt"), barwidth = unit(6, "pt"))
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.12))) +
    panel_text_theme +
    theme(
      axis.text.y = element_text(size = 11, lineheight = 0.95),
      legend.position = "right",
      legend.box = "vertical",
      legend.spacing.y = unit(2, "pt"),
      plot.margin = margin(8, 16, 8, 8)
    )



p.geneset

library(ggpubr)

p.all <- ggarrange(ggarrange(p1, p, nrow=1), ggarrange(p.fpc, p.geneset,nrow=1), nrow=2, ncol=1)

mgam$results |> filter(peak_r_location > 0.7) |> arrange(desc(peak.r)) |> head(n=20)

gene <- "hsp70l"
gene_ix <- match(gene, rownames(mat_r))
if (is.na(gene_ix)) {
  stop("Gene not found in mat_r: ", gene)
}

plot(fit$xyt$r, mat_r[gene_ix,])
plot(fit$xyt$r, mgam$fxs.r[gene_ix,])
plot(fit$xyt$t, mat_r[gene_ix,])
plot(fit$xyt$t, mgam$fxs.t[gene_ix,])
p.hsp.raw <- data.frame(x=fit$xyt$x, y=fit$xyt$y, color=mat_r[gene_ix,]) |>
  ggplot(aes(x=x,y=y,color=color)) + geom_point() +
  scale_color_gradient(low="lightgrey", high="red") + theme_void() + 
  ggtitle(paste0(gene, " expression")) + labs(color="Raw count") +
  guides(color = guide_colorbar(barheight = unit(35, "pt"), barwidth = unit(5, "pt"))) +
  spatial_text_theme

p.hsp.r <- MorphoGAM:::plotGAMestimates(mat_r,
genes=c("hsp70l"),
curve_fit=fit,
mgam=mgam,
type="r"
)
p.hsp.r <- p.hsp.r +
  panel_text_theme +
  theme(
    strip.text = element_text(size = 12),
    axis.title.y = element_text(size = 11, margin = margin(r = 4)),
    plot.margin = margin(8, 8, 8, 8)
  )

p.all <- ggarrange(ggarrange(p1, p, nrow=1, labels=c("a","b"),
                              widths = c(1.05, 1),
                              font.label = list(size = 16, face = "bold")), 
                    ggarrange(p.fpc, p.geneset, labels=c("c","d"), nrow=1,
                              widths = c(1, 1.65),
                              font.label = list(size = 16, face = "bold")),
                    ggarrange(p.hsp.raw, p.hsp.r, nrow=1, labels=c("e","f"),
                              widths = c(1.05, 1),
                              font.label = list(size = 16, face = "bold")),
                    nrow=3,ncol=1, heights=c(1.18, 1.12, 1),
                    align = "hv"
)

ggsave(p.all, filename="../plots/sample_A_figure.png",
       width=13.5, height=12.5, units="in", dpi=300)
