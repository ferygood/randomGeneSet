.libPaths("/data/scratch2/yaochung41/RLib")
source("functions.R")

library(TEKRABber)

load("primateBrainData.RData")

kznf_infer <- read.csv("kznf_bucket.csv")

# get human data
mm_gene <- mmGene
rownames(mm_gene) <- mmGene$geneID
mm_gene <- mm_gene[,-1]

# convert table to TPM
sample_counts <- colSums(mm_gene)

scaling_factor <- sample_counts / 1e6

mm_gene_tpm <- mm_gene
mm_gene_tpm <- t(t(mm_gene_tpm)/ scaling_factor + 1) * 1e6
mm_gene_tpm <- as.data.frame(mm_gene_tpm)

# select KRAB-ZNFs and random gene sets
mm_gene_tpm_name <- ensIDtoGeneName(mm_gene_tpm, species = "mmulatta")

df_kznfs <- mm_gene_tpm_name[
    rownames(mm_gene_tpm_name) %in% kznf_infer$external_gene_name,]

# KRAB-ZNFs against all genes
corrOrthologTE(
    geneInput = df_kznfs,
    teInput = mm_gene_tpm_name,
    numCore = 8,
    fileDir = "./results_allGene_mm",
    fileName = "kznfs_vs_allgene_corr.csv"
)

###########################
# random select 300 genes #
###########################

set.seed(47)
gene_list <- rownames(mm_gene_tpm_name)
selected_genes <- replicate(100, sample(gene_list, size = 300,
                                        replace = FALSE), simplify = FALSE)
for (i in 1:100){

    gene_300 <- selected_genes[[i]]
    df_gene_300 <- mm_gene_tpm_name[
        rownames(mm_gene_tpm_name) %in% selected_genes[[i]], ]

    corrOrthologTE(
        geneInput = df_gene_300,
        teInput = mm_gene_tpm_name,
        numCore = 8,
        fileDir = "./results_allGene_mm",
        fileName = paste0("gene_", i, "_vs_allgene_corr.csv")
    )

}

