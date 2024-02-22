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

# get human TE data
pp_TE <- ppTE

pp_dup <- pp_TE %>%
    group_by(name) %>%
    filter(n()>1)

pp_sum <- pp_dup %>%
    rowwise() %>%
    mutate(sum = sum(c_across(4:100))) %>%
    group_by(name) %>%
    filter(sum == min(sum))

pp_TE <- pp_TE %>%
    anti_join(pp_sum, by=c("name", "family", "class"))

rownames(pp_TE) <- pp_TE$name
pp_TE <- pp_TE[,-c(1,2,3)]

# convert human TE to TPM
te_count <- colSums(pp_TE)
te_scale <- te_count / 1e6
pp_TE_tpm <- t(t(pp_TE)/ te_scale + 1) * 1e6
pp_TE_tpm <- as.data.frame(pp_TE_tpm)

# select KRAB-ZNFs and random gene sets
pp_gene_tpm_name <- ensIDtoGeneName(pp_gene_tpm, species = "ppaniscus")

df_kznfs <- pp_gene_tpm_name[
    rownames(pp_gene_tpm_name) %in% kznf_infer$external_gene_name,]

# KRAB-ZNFs against all genes
corrOrthologTE(
    geneInput = df_kznfs[1:10, ],
    teInput = pp_TE_tpm[1:10, ],
    fileDir = "./results_pp",
    fileName = "kznfs_vs_TE_corr.csv"
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
        geneInput = df_gene_300[1:10, ],
        teInput = pp_TE_tpm[1:10, ],
        fileDir = "./results_pp",
        fileName = paste0("gene_", i, "_vs_TE_corr.csv")
    )

}

