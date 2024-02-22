library(dplyr)

ensIDtoGeneName <- function(df, species="hsapiens"){

    ensIDList <- rownames(df)
    geneName <- gprofiler2::gconvert(
        query=ensIDList,
        organism = species,
        target = "ENSG",
        mthreshold = Inf,
        filter_na = TRUE
    )

    df$ensID <- rownames(df)
    df <- df %>%
        inner_join(geneName[,c("target", "name")], join_by(ensID==target))

    df <- df[!duplicated(df$name), ]

    rownames(df) <- df$name
    end_idx <- ncol(df) - 2
    df <- df[,c(1:end_idx)]

    df

}


