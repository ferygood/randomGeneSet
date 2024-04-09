#.libPaths("/data/scratch2/yaochung41/RLib")
source("functions.R")

library(TEKRABber)
library(dplyr)
library(twice)

load("primateBrainData.RData")
data("hmKZNFs337")

# genes
df_mm_gene <- mmGene[,c(-1)]
rownames(df_mm_gene) <- mmGene$geneID

# transposable elements
mmTE_temp <- mmTE[!duplicated(mmTE$name), ]
mmTEexp <- mmTE_temp %>% select(-c(1,2,3))
rownames(mmTEexp) <- mmTE_temp$name

# genes covert to tpm
sample_counts <- colSums(df_mm_gene)
scaling_factor <- sample_counts / 1e6
df_mm_gene_tpm <- t(t(df_mm_gene)/ scaling_factor + 1) * 1e6
df_mm_gene_tpm <- as.data.frame(df_mm_gene_tpm)
df_mm_kznfs_tpm <- ensIDtoGeneName(df_mm_gene_tpm, species = "mmulatta")
df_mm_kznfs_tpm <- df_mm_kznfs_tpm %>%
    filter(rownames(.) %in% hmKZNFs337$external_gene_name) #292

# tes convert to tpm
te_counts <- colSums(mmTEexp)
te_scale <- te_counts / 1e6
mmTE_tpm <- t(t(mmTEexp)/ scaling_factor + 1) * 1e6
mmTE_tpm <- as.data.frame(mmTE_tpm)


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

mm_c1_corr <- get_corr("cluster1", "Macaca mulatta", df_mm_kznfs_tpm, mmTE_tpm)
mm_c2_corr <- get_corr("cluster2", "Macaca mulatta", df_mm_kznfs_tpm, mmTE_tpm)
mm_c3_corr <- get_corr("cluster3", "Macaca mulatta", df_mm_kznfs_tpm, mmTE_tpm)
mm_c4_corr <- get_corr("cluster4", "Macaca mulatta", df_mm_kznfs_tpm, mmTE_tpm)
mm_c5_corr <- get_corr("cluster5", "Macaca mulatta", df_mm_kznfs_tpm, mmTE_tpm)
mm_c6_corr <- get_corr("cluster6", "Macaca mulatta", df_mm_kznfs_tpm, mmTE_tpm)
mm_c7_corr <- get_corr("cluster7", "Macaca mulatta", df_mm_kznfs_tpm, mmTE_tpm)


mm_tpm_corr <- list(
    "mm_c1_corr" = mm_c1_corr,
    "mm_c2_corr" = mm_c2_corr,
    "mm_c3_corr" = mm_c3_corr,
    "mm_c4_corr" = mm_c4_corr,
    "mm_c5_corr" = mm_c5_corr,
    "mm_c6_corr" = mm_c6_corr,
    "mm_c7_corr" = mm_c7_corr
)

saveRDS(mm_tpm_corr, "mm_tpm_corr.rds")


