#' Create multi-assay experiment object. 
#' 
#' @param study_name study identifier.
#' @return A MultiAssayExperiment object containing clinical and multi-'omics data
#' 
#' @importFrom rlang .data
create_mae <- function(
  # Valid study_name-parameters: tcga, taylor, sun, hieronymus
  study_name = "TCGA",
  # Level of verbosity
  verb = 1,
  ...
){
  if(verb>=1) print(paste("Starting to process:", study_name))
  data_sets <- list.files("data-raw/",pattern=tolower(study_name))
  data_sets <- paste0("data-raw/", data_sets)
  
  # import gene expression matrix from generate.R
  gex_name <- load(grep("gex_",data_sets, value=TRUE))
  gex_object <- get(gex_name) 
  if(verb>=1) print("GEX loaded")
  
  # import gene expression matrix from generate.R
  if (length(grep("cna_",data_sets, value=TRUE))==0) {
    cna_object <- NULL
    if(verb>=1) print("No CNA found")
  } else {
    cna_name <- load(grep("cna_",data_sets, value=TRUE))
    cna_object <- get(cna_name) 
    if(verb>=1) print("CNA loaded")
  }
  
  # import clinical/pheno data from download-clinical.R
  pheno_name <- load(grep("clinical_",data_sets, value=TRUE))
  pheno_object <- get(pheno_name) 
  if(verb>=1) print("Clinical information loaded")
  
  # making sure matrix(s) has samples as columns and genes as rows
  if(length(intersect(row.names(gex_object), pheno_object$sample_name))==0){
    gex_object <- gex_object
  }
  if(!is.null(cna_object) & 
     length(intersect(row.names(cna_object), pheno_object$sample_name))==0){
    cna_object <- cna_object
  }
  
  # need to set pheno_object row names = patient_id for MAE reasons? 
  # with the larger TCGA dataset there is a repeated patient - which makes  this step impossible
  # for the time being - subset to the overlapping patients in gex and clinical? 
  pheno_object <- pheno_object %>% 
    dplyr::filter(.data$sample_name %in% colnames(gex_object))
    
  rownames(pheno_object) <- pheno_object$patient_id

  
  # Generate a MAE-object
  # TODO:
  # Separate here which 'omics are available; e.g. Hieronymus et al. contains only CNA

  if(verb>=1) print("Data reshaping finished")
  
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

  if(verb>=1) print(paste("MAE-object successfully created for", study_name))
  
  return(mae_object)
  
}

# create_tidy <- function(){}
