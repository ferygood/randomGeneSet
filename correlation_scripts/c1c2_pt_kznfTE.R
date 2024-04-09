#.libPaths("/data/scratch2/yaochung41/RLib")
source("functions.R")

library(TEKRABber)
library(dplyr)
library(twice)

load("primateBrainData.RData")
data("hmKZNFs337")

# genes
df_pt_gene <- ptGene[,c(-1)]
rownames(df_pt_gene) <- ptGene$geneID

# transposable elements
ptTEexp <- ptTE %>% select(-c(1,2,3))
rownames(ptTEexp) <- ptTE$name

# genes covert to tpm
sample_counts <- colSums(df_pt_gene)
scaling_factor <- sample_counts / 1e6
df_pt_gene_tpm <- t(t(df_pt_gene)/ scaling_factor + 1) * 1e6
df_pt_gene_tpm <- as.data.frame(df_pt_gene_tpm)
df_pt_kznfs_tpm <- ensIDtoGeneName(df_pt_gene_tpm, species = "ptroglodytes")
df_pt_kznfs_tpm <- df_pt_kznfs_tpm %>%
    filter(rownames(.) %in% hmKZNFs337$external_gene_name) #292

# tes convert to tpm
te_counts <- colSums(ptTEexp)
te_scale <- te_counts / 1e6
ptTE_tpm <- t(t(ptTEexp)/ scaling_factor + 1) * 1e6
ptTE_tpm <- as.data.frame(ptTE_tpm)


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

pt_c1_corr <- get_corr("cluster1", "Pan troglodytes", df_pt_kznfs_tpm, ptTE_tpm)
pt_c2_corr <- get_corr("cluster2", "Pan troglodytes", df_pt_kznfs_tpm, ptTE_tpm)
pt_c3_corr <- get_corr("cluster3", "Pan troglodytes", df_pt_kznfs_tpm, ptTE_tpm)
pt_c4_corr <- get_corr("cluster4", "Pan troglodytes", df_pt_kznfs_tpm, ptTE_tpm)
pt_c5_corr <- get_corr("cluster5", "Pan troglodytes", df_pt_kznfs_tpm, ptTE_tpm)
pt_c6_corr <- get_corr("cluster6", "Pan troglodytes", df_pt_kznfs_tpm, ptTE_tpm)
pt_c7_corr <- get_corr("cluster7", "Pan troglodytes", df_pt_kznfs_tpm, ptTE_tpm)

pt_tpm_corr <- list(
    "pt_c1_corr" = pt_c1_corr,
    "pt_c2_corr" = pt_c2_corr,
    "pt_c3_corr" = pt_c3_corr,
    "pt_c4_corr" = pt_c4_corr,
    "pt_c5_corr" = pt_c5_corr,
    "pt_c6_corr" = pt_c6_corr,
    "pt_c7_corr" = pt_c7_corr
)

saveRDS(pt_tpm_corr, "pt_tpm_corr.rds")


