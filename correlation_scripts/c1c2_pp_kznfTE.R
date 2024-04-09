#.libPaths("/data/scratch2/yaochung41/RLib")
source("functions.R")

library(TEKRABber)
library(dplyr)
library(twice)

load("primateBrainData.RData")
data("hmKZNFs337")

# genes
df_pp_gene <- ppGene[,c(-1)]
rownames(df_pp_gene) <- ppGene$geneID

# transposable elements
pppE_temp <- pppE[!duplicated(pppE$name), ]
pppEexp <- pppE_temp %>% select(-c(1,2,3))
rownames(pppEexp) <- pppE_temp$name

# genes covert to tpm
sample_counts <- colSums(df_pp_gene)
scaling_factor <- sample_counts / 1e6
df_pp_gene_tpm <- t(t(df_pp_gene)/ scaling_factor + 1) * 1e6
df_pp_gene_tpm <- as.data.frame(df_pp_gene_tpm)
df_pp_kznfs_tpm <- ensIDtoGeneName(df_pp_gene_tpm, species = "ppaniscus")
df_pp_kznfs_tpm <- df_pp_kznfs_tpm %>%
    filter(rownames(.) %in% hmKZNFs337$external_gene_name) #292

# tes convert to tpm
te_counts <- colSums(pppEexp)
te_scale <- te_counts / 1e6
pppE_tpm <- t(t(pppEexp)/ scaling_factor + 1) * 1e6
pppE_tpm <- as.data.frame(pppE_tpm)


# select cluster 1 and 2 IDs
cluster_meta <- read.csv("cluster_meta.csv")

get_corr <- function(cluster_num, species, geneTable, teTable){

    cluster_region <- cluster_meta %>%
        filter(cluster==cluster_num) %>%
        select(brain_region) %>%
        unique()

    cluster_id <- metadata %>%
        filter(Organism==species & brain_region %in% cluster_region$brain_region) %>%
        select(1)

    df_gene <- geneTable %>%
        select(cluster_id$Run)

    df_te <- teTable %>%
        select(cluster_id$Run)

    df_corr <- corrOrthologTE(
        geneInput = df_gene,
        teInput = df_te,
        numCore = 3
    )

    df_corr

}

pp_c1_corr <- get_corr("cluster1", "Pan paniscus", df_pp_kznfs_tpm, pppE_tpm)
pp_c2_corr <- get_corr("cluster2", "Pan paniscus", df_pp_kznfs_tpm, pppE_tpm)
pp_c3_corr <- get_corr("cluster3", "Pan paniscus", df_pp_kznfs_tpm, pppE_tpm)
pp_c4_corr <- get_corr("cluster4", "Pan paniscus", df_pp_kznfs_tpm, pppE_tpm)
pp_c5_corr <- get_corr("cluster5", "Pan paniscus", df_pp_kznfs_tpm, pppE_tpm)
pp_c6_corr <- get_corr("cluster6", "Pan paniscus", df_pp_kznfs_tpm, pppE_tpm)
pp_c7_corr <- get_corr("cluster7", "Pan paniscus", df_pp_kznfs_tpm, pppE_tpm)


pp_tpm_corr <- list(
    "pp_c1_corr" = pp_c1_corr,
    "pp_c2_corr" = pp_c2_corr,
    "pp_c3_corr" = pp_c3_corr,
    "pp_c4_corr" = pp_c4_corr,
    "pp_c5_corr" = pp_c5_corr,
    "pp_c6_corr" = pp_c6_corr,
    "pp_c7_corr" = pp_c7_corr
)

saveRDS(pp_tpm_corr, "pp_tpm_corr.rds")



