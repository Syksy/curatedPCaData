#' Create multi-assay experiment object. 
#' 
#' @param study_name study identifier.
#' @return A MulitAssayExperiment object containing clinical and gene expression data 
#' along with any other available data 
#' 
#' @importFrom rlang .data
create_mae <- function(
  study_name = "TCGA"
){
  data_sets <- list.files("data-raw/",pattern=paste0("*_",tolower(study_name),".RData"))
  data_sets <- paste0("data-raw/", data_sets)
  
  # import gene expression matrix from generate.R
  gex_name <- load(grep("data-raw/gex_.*.RData",data_sets, value=TRUE))
  # gex_name <- load(paste0("data-raw/gex_",tolower(study_name),".Rdata"))
  gex_object <- get(gex_name) 
  
  # import gene expression matrix from generate.R
  if (length(grep("data-raw/cna_.*.RData", data_sets, value=TRUE))==0) {
    cna_object <- NULL
  } else {
    cna_name <- load(grep("data-raw/cna_.*.RData", data_sets, value=TRUE))
    cna_object <- get(cna_name) 
  }
  
  # import clinical/pheno data from download-clinical.R
  pheno_name <- load(grep("data-raw/clinical_.*.RData", data_sets, value=TRUE))
  pheno_object <- get(pheno_name) 
  
  # need to set pheno_object row names = patient_id for MAE reasons? 
  # with the larger TCGA dataset there is a repeated patient - which makes  this step impossible
  # for the time being - subset to the overlapping patients in gex and clinical? 
  pheno_object <- pheno_object %>% 
    dplyr::filter(.data$sample_name %in% colnames(gex_object))
    
  rownames(pheno_object) <- pheno_object$patient_id
  
  # Generate a MAE-object
  mae_object <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = if(is.null(cna_object)){
      MultiAssayExperiment::ExperimentList(
        "GEX" = gex_object
      )
    } else{
      MultiAssayExperiment::ExperimentList(
        "GEX" = gex_object,
        "CNA" = cna_object
      )
    },
    colData = pheno_object,
    sampleMap = if(is.null(cna_object)){
      data.frame(assay = "GEX", primary = pheno_object$patient_id,
                 colname = pheno_object$sample_name)
    } else {
      rbind(data.frame(assay = "GEX", primary = pheno_object$patient_id,
                       colname = pheno_object$sample_name),
            data.frame(assay = "CNA", primary = pheno_object$patient_id,
                       colname = pheno_object$sample_name))
    }
    )
  
  return(mae_object)
  
}

# create_tidy <- function(){}
