# Risk Score Functions ----------------------------------------------------

#' Get list of genes to use for Prolaris Risk Score calculation
#'
#' @return list of genes for calculating Prolaris Risk Score
#'
#' @examples
#' prolaris_genes <- getProlarisGenes()
#'
#' @noRd
#' @keywords internal
getProlarisGenes <- function() {
    return(list(
        "FOXM1" = c("FOXM1"),
        "CDC20" = c("CDC20"),
        "CDKN3" = c("CDKN3"),
        "CDC2" = c("CDC2"),
        "KIF11" = c("KIF11"),
        "KIAA0101" = c("KIAA0101"),
        "NUSAP1" = c("NUSAP1"),
        "CENPF" = c("CENPF"),
        "ASPM" = c("ASPM"),
        "BUB1B" = c("BUB1B"),
        "RRM2" = c("RRM2"),
        "DLGAP5" = c("DLGAP5"),
        "BIRC5" = c("BIRC5"),
        "KIF20A" = c("KIF20A"),
        "PLK1" = c("PLK1"),
        "TOP2A" = c("TOP2A"),
        "TK1" = c("TK1"),
        "PBK" = c("PBK"),
        "ASF1B" = c("ASF1B"),
        "C18orf24" = c("C18orf24"),
        "RAD54L" = c("RAD54L"),
        "PTTG1" = c("PTTG1"),
        "CDCA3" = c("CDCA3"),
        "MCM10" = c("MCM10"),
        "PRC1" = c("PRC1"),
        "DTL" = c("DTL"),
        "CEP55" = c("CEP55"),
        "RAD51" = c("RAD51"),
        "CENPM" = c("CENPM"),
        "CDCA8" = c("CDCA8"),
        "ORC6L" = c("ORC6L"),
        "SKA1" = c("SKA1"),
        "ORC6" = c("ORC6"),
        "CDK1" = c("CDK1")
    ))
}

#' Get Genes for Oncotype DX Risk Score calculation
#'
#' @return list of genes used for calculating Oncotype DX Risk Score
#'
#' @examples
#' oncotype_genes <- getOncotypeGenes()
#'
#' @noRd
#' @keywords internal
getOncotypeGenes <- function() {
    return(list(
        "AZGP1" = c("AZGP1"),
        "KLK2" = c("KLK2"),
        "SRD5A2" = c("SRD5A2"),
        "FAM13C" = c("FAM13C"),
        "FLNC" = c("FLNC"),
        "GSN" = c("GSN"),
        "TPM2" = c("TPM2"),
        "GSTM2" = c("GSTM2"),
        "TPX2" = c("TPX2"),
        "BGN" = c("BGN"),
        "COL1A1" = c("COL1A1"),
        "SFRP4" = c("SFRP4")
    ))
}

#' Get "over" Genes for Decipher Risk Score Calculation
#'
#' Genes were pulled from Tabe 2:
#' https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3691249/
#'
#' @return list of genes used for Decipher risk score calculation
#'
#' @examples
#' decipher_over_genes <- getDecipherOverGenes()
#'
#' @noRd
#' @keywords internal
getDecipherOverGenes <- function() {
    return(list(
        "CAMK2N1" = c("CAMK2N1"),
        "EPPK1" = c("EPPK1"),
        "IQGAP3" = c("IQGAP3"),
        "LASP1" = c("LASP1"),
        "NFIB" = c("NFIB"),
        "NUSAP1" = c("NUSAP1"),
        "PBX1" = c("PBX1"),
        "S1PR4" = c("S1PR4"),
        "THBS2" = c("THBS2"),
        "UBE2C" = c("UBE2C"),
        "ZWILCH" = c("ZWILCH")
    ))
}

#' Get "under" Genes for Decipher Risk Score Calculation
#'
#' Genes were pulled from Tabe 2:
#' https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3691249/
#'
#' @return list of genes used for Decipher risk score calculation
#'
#' @examples
#' decipher_over_genes <- getDecipherUnderGenes()
#'
#' @noRd
#' @keywords internal
getDecipherUnderGenes <- function() {
    return(list(
        "ANO7" = c("ANO7"),
        "C6orf10" = c("C6orf10", "TSBP1"),
        "PCDH7" = c("PCDH7"),
        "MYBPC1" = c("MYBPC1"),
        "TSBP" = c("TSBP"),
        "RABGAP1" = c("RABGAP1"),
        "PCAT-32" = c("PCAT-32", "PCAT1"),
        "PCAT-80" = c("PCAT-80", "GLYATL1P4"),
        "TNFRSF19" = c("TNFRSF19")
    ))
}

#' Get a list of genes for AR Risk Score calculation
#'
#' @return list of genes used to calculate AR Risk Score
#'
#' @examples
#' ar_genes <- getARGenes()
#'
#' @noRd
#' @keywords internal
getARGenes <- function() {
    return(list(
        "KLK3" = c("KLK3", "PSA", "APS", "KLK2A1"), # Possibly HK3; ambiguous
        "KLK2" = c("KLK2", "HGK-1", "HGK.1", "KLK2A2"), # Possibly HK2; ambiguous
        "PMEPA1" = c("PMEPA1", "STAG1", "TMEPAI"),
        "ABCC4" = c("ABCC4", "MRP4", "MOATB", "MOAT-B", "MOAT.B"),
        "NKX3-1" = c("NKX3-1", "NKX3.1", "BAPX2", "NKX3A"),
        "C1orf116" = c("C1orf116", "SARG", "FLJ36507", "MGC2742", "MGC4309"),
        "FKBP5" = c("FKBP5", "FKBP51", "FKBP54", "FKBP-51", "FKBP.51", "AIG6", "FKBP-5", "FKBP.5"),
        "ACSL3" = c("ACSL3", "FACL3", "LACS3"),
        "ZBTB10" = c("ZBTB10", "RINZFC", "RINZF"),
        "HERC3" = c("HERC3"),
        "PTGER4" = c("PTGER4", "EP4", "P4R"), # Possibly PTGER2; ambiguous
        "MPHOSPH9" = c("MPHOSPH9", "MPHOS9", "MPP9", "MPP-9", "MPP.9"),
        "EAF2" = c("EAF2", "TRAITS", "BM040", "U19"),
        "MED28" = c("MED28", "EG1"),
        "NNMT" = c("NNMT"),
        "MAF" = c("MAF", "CTRCT21", "AYGRP", "CCA4", "C-MAF", "C.MAF"),
        "GNMT" = c("GNMT"),
        "CENPN" = c("CENPN", "ICEN32", "BM039", "FLJ13607", "FLJ22660"),
        "ELL2" = c("ELL2", "MRCCAT1"),
        "TMPRSS2" = c("TMPRSS2", "PRSS10")
    ))
}

#' Retrieve the gene aliases in BioMart for a gene
#'
#' @param gene single HUGO gene symbol
#'
#' @return vector of gene aliases for input HUGO symbol
#'
#' @examples
#' p53_aliases <- expandAliases("TP53")
#'
#' @noRd
#' @keywords internal
expandAliases <- function(gene) {
    if (length(which(curatedPCaData:::curatedPCaData_genes$hgnc_symbol == gene)) > 0) {
        unique(c(
            gene,
            unlist(
                strsplit(curatedPCaData:::curatedPCaData_genes[
                    which(curatedPCaData:::curatedPCaData_genes$hgnc_symbol == gene), "Aliases"
                ], ";")[[1]]
            )
        ))
    } else {
        gene
    }
}

#' Prepare gene expression matrices for risk score calculations
#'
#' @param GEX gene expression matrix
#' @param gene_list list of genes to prepare
#' @param log_transform whether the gene expression matrix should be log transformed or not
#'
#' @return Modified gene expression matrix
#'
#' @examples
#'
#' GEX <- mae_tcga@ExperimentList$gex.rsem.log
#' gene_list <- lapply(getOncotypeGenes, expandAliases)
#' preppedGEX <- prepGEX(GEX, gene_list, log_transform = FALSE)
#'
#' @noRd
#' @keywords internal
prepGEX <- function(GEX, gene_list, log_transform = TRUE) {
    # check if we get what we are expecting
    if (length(dim(GEX)) != 2) {
        stop("Please use a gene expression matrix with more than 1 dimension")
    }
    if (!inherits(gene_list, "list")) {
        stop("Please enter a gene list")
    }

    # check orientation of gene expression table
    # make column names genes
    tmp_genes <- unlist(gene_list)
    if (TRUE %in% (tmp_genes %in% row.names(GEX))) {
        GEX <- t(GEX)
    }

    # transform log if needed
    if (log_transform) {
        GEX <- log2(GEX + 1)
    }

    GEX <- data.frame(GEX, check.names = F)
    # select main gene, or average alias matches
    GEX_prepped_list <- lapply(seq(gene_list), function(x) {
        if (names(gene_list)[x] %in% colnames(GEX)) {
            return(GEX[names(gene_list)[x]])
        } else {
            dat <- t(data.frame(as.list(rowMeans(GEX[colnames(GEX) %in% gene_list[[x]]], na.rm = T))))
            colnames(dat)[1] <- names(gene_list)[x]
            return(dat)
        }
    })

    # list to df and drop missing genes
    GEX_prepped_df <- cbind.data.frame(GEX_prepped_list)
    GEX_prepped_df[, !apply(is.na(GEX_prepped_df), 2, all)]
}


# Oncotype DX Functions ---------------------------------------------------

cellular_organization_module <- function(GEX_df) {
    val <- (0.163 * (if (is.null(GEX_df$FLNC)) rep(0, nrow(GEX_df)) else GEX_df$FLNC)) +
        (0.504 * (if (is.null(GEX_df$GSN)) rep(0, nrow(GEX_df)) else GEX_df$GSN)) +
        (0.421 * (if (is.null(GEX_df$TPM2)) rep(0, nrow(GEX_df)) else GEX_df$TPM2)) +
        (0.394 * (if (is.null(GEX_df$GSTM2)) rep(0, nrow(GEX_df)) else GEX_df$GSTM2))
    return(val)
}
stromal_module <- function(GEX_df) {
    val <- (0.527 * (if (is.null(GEX_df$BGN)) rep(0, nrow(GEX_df)) else GEX_df$BGN)) +
        (0.457 * (if (is.null(GEX_df$COL1A1)) rep(0, nrow(GEX_df)) else GEX_df$COL1A1)) +
        (0.156 * (if (is.null(GEX_df$SFRP4)) rep(0, nrow(GEX_df)) else GEX_df$SFRP4))
    return(val)
}
androgen_module <- function(GEX_df) {
    val <- (0.634 * (if (is.null(GEX_df$FAM13C)) rep(0, nrow(GEX_df)) else GEX_df$FAM13C)) +
        (1.079 * (if (is.null(GEX_df$KLK2)) rep(0, nrow(GEX_df)) else GEX_df$KLK2)) +
        (0.997 * (if (is.null(GEX_df$SRD5A2)) rep(0, nrow(GEX_df)) else ifelse(GEX_df$SRD5A2 < 5.5, 5.5, GEX_df$SRD5A2))) +
        (0.642 * (if (is.null(GEX_df$AZGP1)) rep(0, nrow(GEX_df)) else GEX_df$AZGP1))
    return(val)
}
proliferation_module <- function(GEX_df) {
    val <- (if (is.null(GEX_df$TPX2)) rep(0, nrow(GEX_df)) else ifelse(GEX_df$TPX2 < 5.0, 5.0, GEX_df$TPX2))
    return(val)
}
#' Oncotype DX Score calculation
#'
#' @param GEX_df gene expression matrix prepped with prepGEX
#'
#' @return vector of scores
#'
#' @examples
#'
#' GEX <- mae_tcga@ExperimentList$gex.rsem.log
#' gene_list <- lapply(getOncotypeGenes, expandAliases)
#' preppedGEX <- prepGEX(GEX, gene_list, log_transform = FALSE)
#' oncotype <- oncotype_score(preppedGEX)
#'
#' @noRd
#' @keywords internal
oncotype_score <- function(GEX_df) {
    risk_score <- 0.735 * stromal_module(GEX_df) -
        0.368 * cellular_organization_module(GEX_df) -
        0.352 * androgen_module(GEX_df) +
        0.095 * proliferation_module(GEX_df)
    names(risk_score) <- row.names(GEX_df)
    return(risk_score)
}

# Prolaris Functions ------------------------------------------------------

#' Prolaris score calculation
#'
#' @param GEX_df gene expression matrix prepped with prepGEX
#'
#' @return vector of scores
#'
#' @examples
#'
#' GEX <- mae_kunderfranco@ExperimentList$gex.logr
#' gene_list <- lapply(getProlarisGenes, expandAliases)
#' preppedGEX <- prepGEX(GEX, gene_list, log_transform = FALSE)
#' prolaris <- prolaris_score(preppedGEX)
#'
#' @noRd
#' @keywords internal
prolaris_score <- function(GEX_df) {
    medians <- apply(GEX_df, 2, stats::median)
    centered_dat <- sweep(GEX_df, 2, medians)
    centered_dat_sq <- centered_dat^2

    score <- apply(centered_dat_sq, 1, mean) %>%
        log2()

    return(score)
}

# Decipher Functions ------------------------------------------------------

#' Decipher score calculation
#'
#' @param GEX_df gene expression matrix prepped with prepGEX
#'
#' @return vector of scores
#'
#' @examples
#'
#' GEX <- mae_tcga@ExperimentList$gex.rsem.log
#' gene_list <- lapply(c(
#'     decipher_genes_over,
#'     decipher_genes_under
#' ), expandAliases)
#' preppedGEX <- prepGEX(GEX, gene_list, log_transform = FALSE)
#' decipher <- decipher_score(preppedGEX)
#'
#' @noRd
#' @keywords internal
decipher_score <- function(GEX_df) {
    over_genes <- unlist(getDecipherOverGenes())
    under_genes <- unlist(getDecipherUnderGenes())
    medians <- apply(GEX_df, 2, median)
    centered_dat <- sweep(GEX_df, 2, medians)
    c1 <- apply(centered_dat[, intersect(over_genes, colnames(centered_dat))], 1, mean)
    c2 <- apply(centered_dat[, intersect(under_genes, colnames(centered_dat))], 1, mean)
    return(c1 - c2)
}

# AR Functions ------------------------------------------------------------

#' AR score calculation
#'
#' @param GEX_df gene expression matrix prepped with prepGEX
#'
#' @return vector of scores
#'
#' @examples
#'
#' GEX <- mae_tcga@ExperimentList$gex.rsem.log
#' gene_list <- lapply(getARGenes, expandAliases)
#' preppedGEX <- prepGEX(GEX, gene_list, log_transform = FALSE)
#' ar <- at_score(preppedGEX)
#'
#' @noRd
#' @keywords internal
ar_score <- function(GEX_df) {
    rowSums(scale(GEX_df))
}
