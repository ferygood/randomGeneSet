.libPaths("/data/scratch2/yaochung41/RLib")
source("functions.R")

library(TEKRABber)

load("primateBrainData.RData")

kznf_infer <- read.csv("kznf_bucket.csv")

# get human data
hm_gene <- hmGene
rownames(hm_gene) <- hmGene$geneID
hm_gene <- hm_gene[,-1]

# convert table to TPM
sample_counts <- colSums(hm_gene)

scaling_factor <- sample_counts / 1e6

hm_gene_tpm <- hm_gene
hm_gene_tpm <- t(t(hm_gene_tpm)/ scaling_factor + 1) * 1e6
hm_gene_tpm <- as.data.frame(hm_gene_tpm)

# get human TE data
hm_TE <- hmTE
rownames(hm_TE) <- hmTE$name
hm_TE <- hm_TE[,-c(1,2,3)]

# convert human TE to TPM
te_count <- colSums(hm_TE)
te_scale <- te_count / 1e6
hm_TE_tpm <- t(t(hm_TE)/ te_scale + 1) * 1e6
hm_TE_tpm <- as.data.frame(hm_TE_tpm)

# select KRAB-ZNFs and random gene sets
hm_gene_tpm_name <- ensIDtoGeneName(hm_gene_tpm)

df_kznfs <- hm_gene_tpm_name[
    rownames(hm_gene_tpm_name) %in% kznf_infer$external_gene_name,]

# KRAB-ZNFs against all genes
corrOrthologTE(
    geneInput = df_kznfs[1:10, ],
    teInput = hm_TE_tpm[1:10, ],
    fileDir = "./results",
    fileName = "kznfs_vs_TE_corr.csv"
)

###########################
# random select 300 genes #
###########################

set.seed(47)
gene_list <- rownames(hm_gene_tpm_name)
selected_genes <- replicate(100, sample(gene_list, size = 300,
                                        replace = FALSE), simplify = FALSE)
for (i in 1:100){

    gene_300 <- selected_genes[[i]]
    df_gene_300 <- hm_gene_tpm_name[
        rownames(hm_gene_tpm_name) %in% selected_genes[[i]], ]

    corrOrthologTE(
        geneInput = df_gene_300[1:10, ],
        teInput = hm_TE_tpm[1:10, ],
        fileDir = "./results",
        fileName = paste0("gene_", i, "_vs_TE_corr.csv")
    )

}

