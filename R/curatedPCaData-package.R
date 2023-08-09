#' curatedPCaData
#'
#' The package curatedPCaData offers a selection of annotated
#' prostate cancer datasets featuring multiple omics, manually
#' curated metadata, and derived downstream variables. The studies
#' are offered as MultiAssayExperiment (MAE) objects via ExperimentHub,
#' and comprise of clinical characteristics tied to gene expression,
#' copy number alteration and somatic mutation data. Further, downstream
#' features computed from these multi-omics data are offered. Multiple
#' vignettes help grasp characteristics of the various studies and provide
#' example exploratory and meta-analysis of leveraging the multiple
#' studies provided here-in.
#'
#' @docType package
#' @name curatedPCaData
#' @keywords internal
#' @import S4Vectors
#' @import MultiAssayExperiment
#' @import RaggedExperiment
#' @import ExperimentHub
#' @importFrom methods slot is
#' @importFrom utils read.table data read.csv
#' @importFrom rlang .data
#' @importFrom stats median
#' @importFrom AnnotationHub query
"_PACKAGE"

# The following block is used by usethis to automatically manage roxygen
# namespace tags. Modify with care!  usethis namespace: start usethis
# namespace: end
NULL
