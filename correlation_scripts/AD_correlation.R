.libPaths("/data/scratch2/yaochung41/RLib")
library(TEKRABber)
library(dplyr)

# load synapse data
load("mayoCorrInput.Rda")

tcx_gene <- mayo$tcxDE$geneCorrInputRef
tcx_te <- mayo$tcxDE$teCorrInputRef

cbe_gene <- mayo$cbeDE$geneCorrInputRef
cbe_te <- mayo$cbeDE$teCorrInputRef

kznfs <- read.csv("kznf_bucket.csv")

# kznfs to TEs
tcx_kznf <- tcx_gene %>%
    filter(rownames(.) %in% kznfs$external_gene_name)

cbe_kznf <- cbe_gene %>%
    filter(rownames(.) %in% kznfs$external_gene_name)

corrOrthologTE(
    geneInput = tcx_kznf,
    teInput = tcx_te,
    numCore = 3,
    fileDir = "./results_tcx",
    fileName = "kznfs_vs_TE_corr.csv"
)

corrOrthologTE(
    geneInput = cbe_kznf,
    teInput = cbe_te,
    numCore = 3,
    fileDir = "./results_cbe",
    fileName = "kznfs_vs_TE_corr.csv"
)

# random genes to TEs
set.seed(11)
gene_list <- rownames(tcx_gene)
select_gene <- replicate(1000, sample(gene_list, size=300,
                                      replace=FALSE), simplify = FALSE)

for (i in 1:2){
    gene_300 <- select_gene[[i]]
    df_tcx <- tcx_gene[rownames(tcx_gene) %in% select_gene[[i]], ]
    df_cbe <- cbe_gene[rownames(cbe_gene) %in% select_gene[[i]], ]

    corrOrthologTE(
        geneInput = df_tcx,
        teInput = tcx_te,
        numCore = 3,
        fileDir = "./results_tcx",
        fileName = paste0("gene_", i, "_vs_TE_corr.csv")
    )

    corrOrthologTE(
        geneInput = df_cbe,
        teInput = cbe_te,
        numCore = 3,
        fileDir = "./results_cbe",
        fileName = paste0("gene_", i, "_vs_TE_corr.csv")
    )
}
