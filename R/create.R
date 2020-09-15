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
  verb = TRUE,
  ...
){
  if(verb == TRUE) print(paste("Starting to process:", study_name))
  data_sets <- list.files("data-raw/",pattern=tolower(study_name))
  if(length(data_sets) == 0){
    stop(paste0("Data for study name ", study_name, " not found; please check spelling"))
  }
  
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
  if(verb == TRUE) print(paste0("Omics: ", paste(omics_names, collapse = ", ")))
  
  ## Construct whole path to data sets
  data_sets <- paste0("data-raw/", data_sets)  
  
  pheno_name <- load(grep("data-raw/clinical_.*.RData", data_sets, value=TRUE))
  pheno_object <- get(pheno_name) 
  if(verb == TRUE) print("Clinical information found")

  # suggest we also pull the sample_name - MAE object doesnt want it here anyways 
  clinical_object <- pheno_object %>%
    dplyr::distinct(.data$patient_id)

  # this line defeats the whole purpose of the MAE
  # rownames(pheno_object) <- make.unique(pheno_object$patient_id)

  # Construct a generalizable sample map
  map <- data.frame(assay = character(0), primary = character(0), colname = character(0))
  for(i in 1:length(omics)){
  	map <- rbind(map, data.frame(assay = names(omics)[[i]], 
  	                             #primary = pheno_object$patient_id, 
  	                             colname = colnames(omics[[i]]))
  	             )
  }
  
  map <- map %>% 
    dplyr::left_join(pheno_object %>% dplyr::select(primary = patient_id, sample_name), 
                     by = c("colname" = "sample_name")) 
  
  # Generate a MAE-object, generalization to various data compositions done above
  if(verb == TRUE) print("Final reformatting")  
  mae_object <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = omics,
    colData = pheno_object,
    sampleMap = map
  )

  if(verb == TRUE) print(paste("MAE-object successfully created for", study_name))
  
  return(mae_object)
}

# create_tidy <- function(){}
