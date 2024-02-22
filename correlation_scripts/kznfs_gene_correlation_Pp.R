.libPaths("/data/scratch2/yaochung41/RLib")
source("functions.R")

library(TEKRABber)

load("primateBrainData.RData")

kznf_infer <- read.csv("kznf_bucket.csv")

# get human data
pp_gene <- ppGene
rownames(pp_gene) <- ppGene$geneID
pp_gene <- pp_gene[,-1]

# convert table to TPM
sample_counts <- colSums(pp_gene)

scaling_factor <- sample_counts / 1e6

pp_gene_tpm <- pp_gene
pp_gene_tpm <- t(t(pp_gene_tpm)/ scaling_factor + 1) * 1e6
pp_gene_tpm <- as.data.frame(pp_gene_tpm)

# select KRAB-ZNFs and random gene sets
pp_gene_tpm_name <- ensIDtoGeneName(pp_gene_tpm, species = "ppaniscus")

df_kznfs <- pp_gene_tpm_name[
    rownames(pp_gene_tpm_name) %in% kznf_infer$external_gene_name,]

# KRAB-ZNFs against all genes
corrOrthologTE(
    geneInput = df_kznfs,
    teInput = pp_gene_tpm_name,
    numCore = 8,
    fileDir = "./results_allGene_pp",
    fileName = "kznfs_vs_allgene_corr.csv"
)

###########################
# random select 300 genes #
###########################

set.seed(47)
gene_list <- rownames(pp_gene_tpm_name)
selected_genes <- replicate(100, sample(gene_list, size = 300,
                                        replace = FALSE), simplify = FALSE)
for (i in 1:100){

    gene_300 <- selected_genes[[i]]
    df_gene_300 <- pp_gene_tpm_name[
        rownames(pp_gene_tpm_name) %in% selected_genes[[i]], ]

    corrOrthologTE(
        geneInput = df_gene_300,
        teInput = pp_gene_tpm_name,
        numCore = 8,
        fileDir = "./results_allGene_pp",
        fileName = paste0("gene_", i, "_vs_allgene_corr.csv")
    )

}

