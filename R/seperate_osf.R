#' Seperates out just the HGNC gene names, in case of duplicate genes, the one with the highest mean across the samples is kept and the others are removed
#'
#' @param osf_data TPM normalized data from the OSF repository
#'

seperate_osf <- function(osf_data){
  
  #osf_data <- data.table::fread("https://osf.io/m5nh6/download")
  
  osf_data <- rio::import("data-raw/TCGA_PRAD_tpm.tsv")
  
  osf_data_2 <- osf_data %>% tidyr::separate(V1,into =c("a","B","C","D","E","F","G","H","i","j","k","l","m","n","o","p"),sep="[|]")

  osf_data_2 <- dplyr::select(osf_data_2,-c("a","B","C","D","E",'i',"G","H","j","k","l","m","n","o","p")) 

  osf_data_2$testMean <- rowMeans(osf_data_2[,2:559], na.rm=TRUE)
  
  osf_data_2 <- osf_data_2 %>% dplyr::group_by(F) %>% dplyr::top_n(1, testMean)
  
  osf_data_2 <- osf_data_2[,-560]
  
  osf_data_2 <- osf_data_2[!duplicated(osf_data_2$F),]
  osf_data_2 <- textshape::column_to_rownames(osf_data_2, loc = 1)
  
  return(osf_data_2)}