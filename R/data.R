#' TCGA MAE-object
#'
#' MultiAssayExperiment object containing GEX and CNA from TCGA 
#' 
#' @format An MAE object spanning 333 men with prostate cancer
#' \describe{
#'   \item{GEX}{matrix with 19985 rows and 333 columns}
#'   \item{CNA}{matrix with 21761 rows and 333 columns}
#' }
#' 
#' @source \url{https://www.cbioportal.org/study/summary?id=prad_tcga_pub}
#'
"mae_tcga"

#' Sun et al. MAE-object
#'
#' MultiAssayExperiment object containing GEX and CNA from Sun et al.
#' 
#' #' @format An MAE object spanning 79 men 
#' \describe{
#'   \item{GEX}{matrix with 12057 rows and 79 columns}
#' }
#' 
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE25136}
#'
"mae_sun"


