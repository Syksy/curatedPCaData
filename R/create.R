#' Create multi-assay experiment object. 
#' 
#' @param study_name study identifier.
#' @return A MulitAssayExperiment object containing clinical and gene expression data 
#' along with any other available data 
create_mae <- function(
  study_name = "TCGA"
){
  data_sets <- list.files("data/",pattern=paste0(tolower(study_name),"_"))
  data_sets <- paste0("data/", data_sets)
  
  # import gene expression matrix from generate.R
  gex_name <- load(grep("*_gex.rda",data_sets, value=TRUE))
  gex_object <- get(gex_name) 
  
  # import gene expression matrix from generate.R
  if (isEmpty(grep("*_cna.rda",data_sets, value=TRUE))) {
    cna_object <- NULL
  } else {
    cna_name <- load(grep("*_cna.rda",data_sets, value=TRUE))
    cna_object <- get(cna_name) 
  }
  
  # import clinical/pheno data from download-clinical.R
  pheno_name <- load(grep("*_clinical.rda",data_sets, value=TRUE))
  pheno_object <- get(pheno_name) 
  
  # making sure matrix(s) has samples as columns and genes as rows
  if(length(intersect(colnames(gex_object), pheno_object$sample_name))==0){
    gex_object <- t(gex_object)
  }
  if(!is.null(cna_object) & 
     length(intersect(colnames(cna_object), pheno_object$sample_name))==0){
    cna_object <- t(cna_object)
  }
  
  # need to set pheno_object row names = patient_id for MAE reasons? 
  # with the larger TCGA dataset there is a repeated patient - which makes  this step impossible
  # for the time being - subset to the overlapping patients in gex and clinical? 
  pheno_object <- pheno_object %>% 
    dplyr::filter(sample_name %in% colnames(gex_object))
    
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
