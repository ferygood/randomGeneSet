#.libPaths("/data/scratch2/yaochung41/RLib")
#source("functions.R")

library(TEKRABber)
library(dplyr)

load("primateBrainData.RData")

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

# tes convert to tpm
te_counts <- colSums(ptTEexp)
te_scale <- te_counts / 1e6
ptTE_tpm <- t(t(ptTEexp)/ scaling_factor + 1) * 1e6
ptTE_tpm <- as.data.frame(ptTE_tpm)


# select cluster 1 and 2 IDs
cluster_meta <- read.csv("cluster_meta.csv")

c1_region <- cluster_meta %>%
    filter(cluster=="cluster1") %>%
    select(brain_region) %>%
    unique()

c2_region <- cluster_meta %>%
    filter(cluster=="cluster2") %>%
    select(brain_region) %>%
    unique()

c1_id <- metadata %>%
    filter(Organism=="Pan troglodytes" & brain_region %in% c1_region$brain_region) %>%
    select(1)

c2_id <- metadata %>%
    filter(Organism=="Pan troglodytes" & brain_region %in% c2_region$brain_region) %>%
    select(1)

# 1000 iteration
set.seed(47)
gene_list <- rownames(df_pt_gene_tpm)
selected_genes <- replicate(1000,
                            sample(gene_list, size=337, replace = FALSE),
                            simplify = FALSE)

df_temp_c1 <- data.frame(
    neg_count = as.numeric(),
    pos_count = as.numeric()
)

df_temp_c2 <- data.frame(
    neg_count = as.numeric(),
    pos_count = as.numeric()
)

for (i in 1:1000){
    gene <- df_pt_gene_tpm
    te <- ptTE_tpm
    gene_set <- selected_genes[[i]]

    gene_c1 <- gene %>%
        select(c1_id$Run) %>%
        filter(rownames(.) %in% gene_set)

    te_c1 <- te %>%
        select(c1_id$Run)

    gene_c2 <- gene %>%
        select(c2_id$Run) %>%
        filter(rownames(.) %in% gene_set)

    te_c2 <- te %>%
        select(c2_id$Run)

    df_c1 <- corrOrthologTE(
        geneInput = gene_c1,
        teInput = te_c1,
        numCore = 3
    )

    df_c2 <- corrOrthologTE(
        geneInput = gene_c2,
        teInput = te_c2,
        numCore = 3
    )

    df_c1_sig <- df_c1 %>% filter(padj<0.01)
    c1_negCount <- df_c1_sig %>% filter(coef<0) %>% nrow()
    c1_posCount <- df_c1_sig %>% filter(coef>0) %>% nrow()
    #file_c1 <- paste0("results_c1_pt/gene_", i, "_vs_TE_corr.csv")
    #write.csv(df_c1_sig, file=file_c1, row.names = F)

    df_c2_sig <- df_c2 %>% filter(padj<0.01)
    c2_negCount <- df_c2_sig %>% filter(coef<0) %>% nrow()
    c2_posCount <- df_c2_sig %>% filter(coef>0) %>% nrow()
    #file_c2 <- paste0("results_c2_pt/gene_", i, "_vs_TE_corr.csv")
    #write.csv(df_c2_sig, file=file_c2, row.names = F)

    c1_temp <- data.frame(
        neg_count = c1_negCount,
        pos_count = c1_posCount
    )

    df_temp_c1 <- rbind(df_temp_c1, c1_temp)

    c2_temp <- data.frame(
        neg_count = c2_negCount,
        pos_count = c2_posCount
    )

    df_temp_c2 <- rbind(df_temp_c2, c2_temp)
}

write.csv(df_temp_c1, file="pt_c1_result.csv", row.names=F)
write.csv(df_temp_c2, file="pt_c2_result.csv", row.names=F)
