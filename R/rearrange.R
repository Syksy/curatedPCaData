#' Transposes the OSF gene expression matrix containing HGNC genes and converts it into a dataframe
#'
#' @param osf_data_2 the modified OSF matrix that contains the HGNC gene names
#'

rearrange <- function(osf_data_2){
  osf_data_t <- t(osf_data_2)
  osf_data_t <- data.frame(osf_data_t)
  
  osf_data_t <- tibble::rownames_to_column(osf_data_t, "Samples")
  #colnames(osf_data_t) <- osf_data_t[1,]
  #osf_data_t <- osf_data_t[-1, ] 
  
  return (osf_data_t)}
