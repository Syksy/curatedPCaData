#' @title Genes used when fetching data for the curatedPCaData-package
#' @description Genes, e.g. Hugo symbols, used when constructing the curatedPCaData package. This data frame is essential for example when calling the genes in cBioPortal-wrapper. This data.frame contains excess information, e.g. one may wish to limit to using unique/sorted hugo symbols (field 'hgnc_symbol').
#' @format A data frame with 268931 rows and 8 variables:
#' \describe{
#'   \item{\code{ensembl_gene_id}}{character EMSEMBL gene ids}
#'   \item{\code{ensembl_transcript_id}}{character ENSEMBL transcript ids}
#'   \item{\code{hgnc_symbol}}{character Hugo symbols}
#'   \item{\code{refseq_mrna}}{character Refseq NM_##### ids}
#'   \item{\code{chromosome_name}}{character Chromosome location}
#'   \item{\code{start_position}}{integer Starting bp location}
#'   \item{\code{end_position}}{integer Ending bp location}
#'   \item{\code{description}}{character Character string description for the gene.} 
#'}
#' @details This data.frame contains information for the package's internal use, e.g. when calling cBioPortal's wrapper with Hugo symbols. Further, attr-field 'date' shows the exact date of retrieving this data.frame using biomaRt.
"curatedPCaData_genes"
