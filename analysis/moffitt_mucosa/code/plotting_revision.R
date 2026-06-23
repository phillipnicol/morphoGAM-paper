
library(MorphoGAM)
library(RColorBrewer)

setwd(here::here("analysis/moffitt_mucosa/code"))

Y <- read.csv("../../data/GL2_distal_colon_cell_by_gene_raw.csv")

rownames(Y) <- Y[,1]

Y <- Y[,-1]

library(Matrix)

Y <- as.matrix(Y)

Y <- t(Y) #Transpose to get genes x cells

meta <- read.csv("../../data/GL2_distal_colon_cell_type_and_locations_2023.08.11.csv")

library(tidyverse)
meta <- as.data.frame(meta)
meta.sub <- meta |> filter(slice_full_name == "20220518_WT_dcol_slice_3") |>
  filter(spatial_neighborhood_v1 == "Mucosa") |>
  filter(leiden_combined_v2 == "Enterocyte")

Y.sub <- Y[,meta.sub$X]

Y.sub <- Y.sub[rowSums(Y.sub) >= 10,]

old.rownames <- rownames(Y.sub)

#Test expression for these genes
library(SPARK)
locus <- as.matrix(meta.sub[,c("x","y")])
res <- SPARK::sparkx(count_in = Y.sub,
                     locus_in = locus)

spark <- res$res_mtest |> as.data.frame() |>
  arrange(adjustedPval)

#Load gaston
gaston_rank <- read.csv("../../python_models/data/gaston_gene_rankings.csv")
#ranked.list <- gaston_rank$abs_slope
gaston.list <- -log10(gaston_rank$p_value)
gaston.list[is.na(gaston.list)] <- 0
names(gaston.list) <- gaston_rank$gene

library(dplyr)

spateo_res <- read.csv("../../python_models/data/spateo_colon_svg_rankings.csv")

spateo_rank <- data.frame(gene = rownames(Y.sub)) %>%
  left_join(spateo_res, by = "gene") %>%
  mutate(
    pval = case_when(
      is.na(pval) ~ 1,
      pval == 0 ~ 10^(-100),
      TRUE ~ pval
    ),
    missing_from_spateo = is.na(pval) | !(gene %in% spateo_res$gene)
  ) %>%
  arrange(pval)
spateo.list <- -log10(spateo_rank$pval)
names(spateo.list) <- spateo_rank$gene




library(ggrepel)
#results_df <- readRDS("../data/loop_svg.RDS")
nnsvg_res <- readRDS("../data/nnSVG_results.RDS")

#angle <- results_df$angle
#radial.response <- results_df$radial.response
#nnsvg_rank <- nnsvg_res$rank
#p.vals <- spark$adjustedPval

load("../data/mucosa_mgam.RData")
load("../data/curve.RData")

rank_label <- function(gene, spark_rank, nnsvg_rank, spateo_rank, gaston_rank) {
  paste0(
    gene, "\n",
    "Spark rank = ", spark_rank, "   nnSVG rank = ", nnsvg_rank, "\n",
    "Spateo rank = ", spateo_rank, "   GASTON rank = ", gaston_rank
  )
}

gam_panel_theme <- theme(
  strip.text = element_text(size = 8, lineheight = 0.95, margin = margin(3, 2, 3, 2)),
  axis.title = element_text(size = 8),
  axis.text = element_text(size = 7),
  panel.spacing = unit(0.8, "lines")
)

rownames(Y.sub) <- old.rownames
top5 <- order(mgam$results$peak.t, decreasing=TRUE)[1:6]
#top5 <- result[1:5]
for(i in top5) {
  spark_rank <- which(rownames(spark) == rownames(Y.sub)[i])
  nnsvg_rank <- which(rownames(nnsvg_res) == rownames(Y.sub)[i])
  spateo_rank <- which(names(spateo.list) == rownames(Y.sub)[i])
  gaston_rank <- which(names(gaston.list) == rownames(Y.sub)[i])
  rownames(Y.sub)[i] <- rank_label(
    old.rownames[i],
    spark_rank,
    nnsvg_rank,
    spateo_rank,
    gaston_rank
  )
}


p.peak <- plotGAMestimates(Y.sub,
                           genes=order(mgam$results$peak.t,decreasing = TRUE)[1:6],
                           curve_fit = fit,
                           mgam_object = mgam,
                           nrow=2) +
  gam_panel_theme +
  scale_y_sqrt()

ggsave(p.peak, filename="../plots/mouse_mucosa_peak_revise.png")

rownames(Y.sub) <- old.rownames
top5 <- order(mgam$results$range.r, decreasing=TRUE)[1:6]
#top5 <- result[1:5]
for(i in top5) {
  spark_rank <- which(rownames(spark) == rownames(Y.sub)[i])
  nnsvg_rank <- which(rownames(nnsvg_res) == rownames(Y.sub)[i])
  spateo_rank <- which(names(spateo.list) == rownames(Y.sub)[i])
  gaston_rank <- which(names(gaston.list) == rownames(Y.sub)[i])
  rownames(Y.sub)[i] <- rank_label(
    old.rownames[i],
    spark_rank,
    nnsvg_rank,
    spateo_rank,
    gaston_rank
  )
}

p.range <- plotGAMestimates(Y.sub,
                            genes=order(mgam$results$range.r,decreasing = TRUE)[1:6],
                            curve_fit = fit,
                            mgam_object = mgam,
                            type="r",
                            nrow=2) +
  gam_panel_theme +
  xlim(0.3,0.7) +
  scale_y_sqrt()

fpc1 <- plotFPCloading(mgam_object=mgam,
                       curve.fit=fit,
                       L=1,num_genes=5)

fpc2 <- plotFPCloading(mgam_object=mgam,
                       curve.fit=fit,
                       L=2,num_genes=4)




df.full <- readRDS("../data/gsea3_results.RDS")

df.full$Method <- factor(
  df.full$Method,
  levels = c("MorphoGAM", "Spark-X", "nnSVG", "GASTON", "Spateo")
)


library(dplyr)
library(ggplot2)
library(stringr)

df_plot <- df.full %>%
  mutate(
    Interferon = str_detect(
      y,
      regex("interferon|\\bIFN\\b|type\\s*[I1]\\s*interferon|type\\s*II\\s*interferon",
            ignore_case = TRUE)
    )
  )

p.sub <- ggplot(df_plot, aes(x = Method, y = x)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90", color = "grey40") +
  geom_jitter(
    data = dplyr::filter(df_plot, !Interferon),
    width = 0.18, height = 0,
    color = "grey50", size = 1.8, alpha = 0.6
  ) +
  geom_jitter(
    data = dplyr::filter(df_plot, Interferon),
    width = 0.18, height = 0,
    color = "red", size = 2.5
  ) +
  labs(
    x = NULL,
    y = "GSEA test statistic (Scaled)",
    title = "",
    subtitle = "Red points = Interferon-related gene sets; Grey = Other gene sets"
  ) +
  theme_bw()








library(ggpubr)
p <- ggarrange(p.peak, p.range, nrow=2, labels=c("a","b"))

p2 <- ggarrange(p, p.sub, ggarrange(fpc1,fpc2, nrow=1), nrow=3,
                heights=c(2,1,1), labels=c("", "c", "d"))


ggsave(p2, filename="../plots/mucosa_updated_gsea_and_pcs_spateo.png",
       width= 4.78*1.5,
       height= 8.18*1.5,
       units="in")
