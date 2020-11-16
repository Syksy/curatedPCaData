#' Takes in modified OSF matrix with HGNC gene names and TCGA sample IDs, converts it into a dataframe and formats the TCGA sample names to match it to the format available in the package
#'
#' @param dropped_results1 
#'

format <- function(dropped_results1){

  osf <- data.frame(dropped_results1)
  #data.table::setDT(osf, keep.rownames = "V1")
  osf$AliquotBarcode <-gsub("-01A.*","",osf$AliquotBarcode)
  osf[2:554,2:58685] <- as.numeric(unlist(osf[2:554,2:58685]))
  return(osf)}