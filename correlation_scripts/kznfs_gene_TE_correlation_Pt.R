.libPaths("/data/scratch2/yaochung41/RLib")
source("functions.R")

library(TEKRABber)

load("primateBrainData.RData")

kznf_infer <- read.csv("kznf_bucket.csv")

# get human data
pt_gene <- ptGene
rownames(pt_gene) <- ptGene$geneID
pt_gene <- pt_gene[,-1]

# convert table to TPM
sample_counts <- colSums(pt_gene)

scaling_factor <- sample_counts / 1e6

pt_gene_tpm <- pt_gene
pt_gene_tpm <- t(t(pt_gene_tpm)/ scaling_factor + 1) * 1e6
pt_gene_tpm <- as.data.frame(pt_gene_tpm)

# get human TE data
pt_TE <- ptTE
rownames(pt_TE) <- ptTE$name
pt_TE <- pt_TE[,-c(1,2,3)]

# convert human TE to TPM
te_count <- colSums(pt_TE)
te_scale <- te_count / 1e6
pt_TE_tpm <- t(t(pt_TE)/ te_scale + 1) * 1e6
pt_TE_tpm <- as.data.frame(pt_TE_tpm)

# select KRAB-ZNFs and random gene sets
pt_gene_tpm_name <- ensIDtoGeneName(pt_gene_tpm, species = "ptroglodytes")

df_kznfs <- pt_gene_tpm_name[
    rownames(pt_gene_tpm_name) %in% kznf_infer$external_gene_name,]

# KRAB-ZNFs against all genes
corrOrthologTE(
    geneInput = df_kznfs[1:10, ],
    teInput = pt_TE_tpm[1:10, ],
    fileDir = "./results_pt",
    fileName = "kznfs_vs_TE_corr.csv"
)

###########################
# random select 300 genes #
###########################

set.seed(47)
gene_list <- rownames(pt_gene_tpm_name)
selected_genes <- replicate(100, sample(gene_list, size = 300,
                                        replace = FALSE), simplify = FALSE)
for (i in 1:100){

    gene_300 <- selected_genes[[i]]
    df_gene_300 <- pt_gene_tpm_name[
        rownames(pt_gene_tpm_name) %in% selected_genes[[i]], ]

    corrOrthologTE(
        geneInput = df_gene_300[1:10, ],
        teInput = pt_TE_tpm[1:10, ],
        fileDir = "./results_pt",
        fileName = paste0("gene_", i, "_vs_TE_corr.csv")
    )

}

