setwd(here::here("analysis/moffitt_mucosa/code"))

library(tidyverse)

library(piano)

#gene_sets <- loadGSC(file="../data/m5.go.bp.v2025.1.Mm.symbols.gmt",
#                     type="gmt")


gene_sets <- loadGSC(file="../data/mh.all.v2025.1.Mm.symbols.gmt",
                     type="gmt")

# View the top enriched Hallmark pathways
#print(head(summary_table_sorted, 10))

nnsvg <- readRDS("../data/nnSVG_results.RDS")

#nnsvg <- nnsvg[1:100,]
ranked.list <- nnsvg$LR_stat
names(ranked.list) <- rownames(nnsvg)
#names(ranked.list) <- toupper(rownames(nnsvg))


gsa_result <- runGSA(
  geneLevelStats = ranked.list,
  geneSetStat = "mean",       # or "mean", "median", "fisher", etc.
  gsc = gene_sets,
  nPerm = 1000                # Increase for stability
)

print(GSAsummaryTable(gsa_result))

summary_table <- GSAsummaryTable(gsa_result)

summary_table_sorted <- summary_table[order(summary_table$`Stat (non-dir.)`,decreasing = TRUE), ]

# View the top enriched Hallmark pathways
print(head(summary_table_sorted, 10))

summary_table <- GSAsummaryTable(gsa_result)

gsea.nnsvg <- summary_table[order(summary_table$`Stat (non-dir.)`,decreasing = TRUE), ]

# View the top enriched Hallmark pathways
#print(head(summary_table_sorted, 10))

df.nnsvg <- data.frame(y=gsea.nnsvg$Name,
                       x=gsea.nnsvg$`Stat (non-dir.)` |> scale(),
                       Method="nnSVG")



spark <- readRDS("../data/spark_results.RDS")



#spark <- spark[1:100,]
ranked.list <- -log10(spark$combinedPval)
#names(ranked.list) <- toupper(rownames(spark))
names(ranked.list) <- rownames(spark)

gsa_result <- runGSA(
  geneLevelStats = ranked.list,
  geneSetStat = "mean",       # or "mean", "median", "fisher", etc.
  gsc = gene_sets,
  nPerm = 1000                # Increase for stability
)

print(GSAsummaryTable(gsa_result))

summary_table <- GSAsummaryTable(gsa_result)

summary_table_sorted <- summary_table[order(summary_table$`Stat (non-dir.)`,decreasing = TRUE), ]

# View the top enriched Hallmark pathways
print(head(summary_table_sorted, 10))

summary_table <- GSAsummaryTable(gsa_result)

gsea.spark <- summary_table[order(summary_table$`Stat (non-dir.)`,decreasing = TRUE), ]

# View the top enriched Hallmark pathways
#print(head(summary_table_sorted, 10))


df.spark <- data.frame(y=gsea.spark$Name,
                       x=gsea.spark$`Stat (non-dir.)` |> scale(),
                       Method="Spark-X")



load("../data/mucosa_mgam.RData")

mgam.df <- mgam$results |> arrange(desc(peak.t))

#mgam.df <- mgam.df[1:100,]

ranked.list <- mgam.df$peak.t
#names(ranked.list) <- toupper(rownames(mgam.df))
names(ranked.list) <- rownames(mgam.df)


gsa_result <- runGSA(
  geneLevelStats = ranked.list,
  geneSetStat = "mean",       # or "mean", "median", "fisher", etc.
  gsc = gene_sets,
  nPerm = 1000                # Increase for stability
)

summary_table <- GSAsummaryTable(gsa_result)

gsea.morphogam <- summary_table[order(summary_table$`Stat (non-dir.)`,decreasing = TRUE), ]


df.mgam <- data.frame(y=gsea.morphogam$Name,
                      x=gsea.morphogam$`Stat (non-dir.)` |> scale(),
                      Method="MorphoGAM")


## Add GASTON
gaston_rank <- read.csv("../../python_models/data/gaston_gene_rankings.csv")
#ranked.list <- gaston_rank$abs_slope
ranked.list <- -log10(gaston_rank$p_value)
ranked.list[is.na(ranked.list)] <- 0
names(ranked.list) <- gaston_rank$gene



gsa_result <- runGSA(
  geneLevelStats = ranked.list,
  geneSetStat = "mean",       # or "mean", "median", "fisher", etc.
  gsc = gene_sets,
  nPerm = 1000                # Increase for stability
)

summary_table <- GSAsummaryTable(gsa_result)

gsea.gaston <- summary_table[order(summary_table$`Stat (non-dir.)`,decreasing = TRUE), ]

df.gaston <- data.frame(y=gsea.gaston$Name,
                        x=gsea.gaston$`Stat (non-dir.)` |> scale(),
                        Method="GASTON")


## Add Spateo 
spateo_rank <- read.csv("../../python_models/data/spateo_colon_svg_rankings.csv")
#Replace 0 pvalues with 10^{-100}
spateo_rank <- spateo_rank %>%
  mutate(
    pval = ifelse(pval == 0, 10^(-100), pval)
  )
ranked.list <- -log10(spateo_rank$pval)
names(ranked.list) <- spateo_rank$gene



gsa_result <- runGSA(
  geneLevelStats = ranked.list,
  geneSetStat = "mean",       # or "mean", "median", "fisher", etc.
  gsc = gene_sets,
  nPerm = 1000                # Increase for stability
)

summary_table <- GSAsummaryTable(gsa_result)

gsea.spateo <- summary_table[order(summary_table$`Stat (non-dir.)`,decreasing = TRUE), ]

df.spateo <- data.frame(y=gsea.spateo$Name,
                        x=gsea.spateo$`Stat (non-dir.)` |> scale(),
                        Method="Spateo")



df.full <- rbind(df.spark, df.mgam, df.nnsvg, df.gaston, df.spateo)

saveRDS(df.full, "../data/gsea3_results.RDS")


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

p <- ggplot(df_plot, aes(x = Method, y = x)) +
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


ggsave(p,filename="../plots/gsea3_results.png")
