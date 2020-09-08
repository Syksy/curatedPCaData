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
  if(length(data_sets)==0){
  	stop(paste0("Data for study name '", study_name, "' not found; please check spelling and 'data-raw'-subfolder."))
  }
  # import gene expression matrix from generate.R
  #gex_name <- load(grep("gex_",data_sets, value=TRUE))
  #gex_object <- get(gex_name) 
  #if(verb>=1) print("GEX loaded")
  #
  # import gene expression matrix from generate.R
  #if (length(grep("cna_",data_sets, value=TRUE))==0) {
  #  cna_object <- NULL
  #  if(verb>=1) print("No CNA found")
  #} else {
  #  cna_name <- load(grep("cna_",data_sets, value=TRUE))
  #  cna_object <- get(cna_name) 
  #  if(verb>=1) print("CNA loaded")
  #}
  
  ## 'omics vary from study to study, read in as many as possible
  omics_sets <- strsplit(data_sets, "_")
  # Omit clinical information character string split from the 'omics portion
  omics_sets <- omics_sets[-which(unlist(lapply(omics_sets, FUN=function(z) { "clinical" %in% z })))]
  # 'omics names is concatenated from all strsplit prior to last element (which is presumably "_STUDYNAME.{rda,RData}")
  omics_names <- unlist(lapply(omics_sets, FUN=function(z) { paste(z[-length(z)], collapse="_") }))
  # Reconstruct paths to omics sets
  omics_sets <- paste0("data-raw/", unlist(lapply(omics_sets, FUN=function(z) { paste(z, collapse="_") })))
  # Construct a list of 'get' objects for the various omics, load them into environment then construct a list
  omics <- list()
  for(f in omics_sets){
  	omic <- load(f)
  	omics[[length(omics)+1]] <- get(omic)
  }
  # Give correct name for each omics list member
  names(omics) <- omics_names
  if(verb>=2) print(paste(c("Omics: ", omics_names), collapse=", "))

  ## Construct whole path to data sets
  data_sets <- paste0("data-raw/", data_sets)  
  ## Clinical variables are required for each study
  # import clinical/pheno data from download-clinical.R
  pheno_name <- load(grep("clinical_",data_sets, value=TRUE))
  pheno_object <- get(pheno_name) 
  if(verb>=1) print("Clinical information loaded")
  
  # making sure matrix(s) has samples as columns and genes as rows
  #if(length(intersect(row.names(gex_object), pheno_object$sample_name))==0){
  #  gex_object <- gex_object
  #}
  #if(!is.null(cna_object) & 
  #   length(intersect(row.names(cna_object), pheno_object$sample_name))==0){
  #  cna_object <- cna_object
  #}
  
  # need to set pheno_object row names = patient_id for MAE reasons? 
  # with the larger TCGA dataset there is a repeated patient - which makes  this step impossible
  # for the time being - subset to the overlapping patients in gex and clinical? 
  #pheno_object <- pheno_object %>% 
  #  dplyr::filter(.data$sample_name %in% colnames(gex_object))
   
  rownames(pheno_object) <- make.unique(pheno_object$patient_id)

  # Construct a generalizable sample map
  map <- data.frame(assay = character(0), primary = character(0), colname = character(0))
  for(i in 1:length(omics)){
  	map <- rbind(map, data.frame(assay = names(omics)[i], primary = pheno_object$patient_id, colname = pheno_object$sample_name))
  }
  
  # Generate a MAE-object, generalization to various data compositions done above
  if(verb>=1) print("Data reshaping finished")  
  mae_object <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = omics,
    colData = pheno_object,
    sampleMap = map
  )

  if(verb>=1) print(paste("MAE-object successfully created for", study_name))
  
  return(mae_object)
}

# create_tidy <- function(){}
