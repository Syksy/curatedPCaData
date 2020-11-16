#' Compares the OSG GEX matrix to a file which contains the TCGA PRAD sample IDs and maps it to the CGHubAnalysisID in the matrix
#'
#' @param PRAD_osf 
#'

keep_only_PRAD <- function(PRAD_osf){
  
  #tcga_map <- data.table::fread("https://osf.io/7qpsg/download")
  
  tcga_map <- rio::import("data-raw/TCGA_ID_MAP.csv")
  tcga_map <- tcga_map %>% dplyr::filter(Disease == "PRAD")
  
  PRAD_osf <-merge(osf_data_t,tcga_map,by.x = 'Samples',by.y = 'CGHubAnalysisID')
  
  dropped_results <- dplyr::select(PRAD_osf,-c('Aliquot_id','Disease','Samples'))
  dropped_results1 <- dropped_results[,c(58685,1:58684)]

  
  
  dropped_results1 <- dropped_results1[!duplicated(dropped_results1$AliquotBarcode),]
  

  return(dropped_results1)}