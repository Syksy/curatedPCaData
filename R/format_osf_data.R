#' Does pre-formatting of OSF dataset so that the genes are in HGNC format and sample names are TCGA IDs
#'
#' @param osf_data 
#'

format_osf_data <- function(osf_data){
  osf_data <- rio::import("../data-raw/TCGA_PRAD_tpm.tsv")
  
  osf_data_2 <- osf_data %>% tidyr::separate(V1,into =c("a","B","C","D","E","F","G","H","i","j","k","l","m","n","o","p"),sep="[|]")
  
  osf_data_2 <- dplyr::select(osf_data_2,-c("a","B","C","D","E",'i',"G","H","j","k","l","m","n","o","p")) 
  
  osf_data_2$testMean <- rowMeans(osf_data_2[,2:559], na.rm=TRUE)
  
  osf_data_2 <- osf_data_2 %>% dplyr::group_by(F) %>% dplyr::top_n(1, testMean)
  
  osf_data_2 <- osf_data_2[,-560]
  
  osf_data_2 <- osf_data_2[!duplicated(osf_data_2$F),]
  osf_data_2 <- textshape::column_to_rownames(osf_data_2, loc = 1)
  
  osf_data_t <- t(osf_data_2)
  osf_data_t <- data.frame(osf_data_t)
  
  osf_data_t <- tibble::rownames_to_column(osf_data_t, "Samples")
  
  tcga_map <- rio::import("../data-raw/TCGA_ID_MAP.csv")
  tcga_map <- tcga_map %>% dplyr::filter(Disease == "PRAD")
  
  PRAD_osf <-merge(osf_data_t,tcga_map,by.x = 'Samples',by.y = 'CGHubAnalysisID')
  
  dropped_results <- dplyr::select(PRAD_osf,-c('Aliquot_id','Disease','Samples'))
  dropped_results1 <- dropped_results[,c(58685,1:58684)]
  
  
  dropped_results1 <- dropped_results1[!duplicated(dropped_results1$AliquotBarcode),]
  
  osf <- data.frame(dropped_results1)
  osf$AliquotBarcode <-gsub("-01A.*","",osf$AliquotBarcode)
  osf[2:554,2:58685] <- as.numeric(unlist(osf[2:554,2:58685]))
  return(osf)}