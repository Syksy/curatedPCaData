#' Create multi-assay experiment object. 
#' 
#' @param study_name study identifier.
#' @param verb print output as MAE is being created 
#' @param ... additional arguments 
#' @return An MAE object containing clinical and multi-'omics data
#' 
#' @importFrom rlang .data
#'
create_mae <- function(
  # Valid study_name-parameters: tcga, taylor, sun, hieronymus, friedrich, barbieri, ren, chandran,igc,wang,kim,abida
 #study_name = c("TCGA", "Taylor", "Sun", "Hieronymus", "Barbieri", "Ren","kim","abida","igc","wang"),
  study_name,
  # Level of verbosity
  verb = TRUE,
  ...
){
  if(verb == TRUE) print(paste("Starting to process:", study_name))
  data_sets <- list.files("data-raw/",pattern=tolower(study_name))
  if(length(data_sets) == 0){
    stop(paste0("Data for study name ", study_name, " not found; please check spelling"))
  }
  
  ## Omit any files that do not contain file type '.RData' 
  if(verb) print(paste("Found data prior to filtering:", paste(data_sets, collapse=", ")))
  data_sets <- data_sets[intersect(1:length(data_sets), grep('.RData', data_sets))]
  if(verb) print(paste("Found data post filtering:", paste(data_sets, collapse=", ")))

  ## 'omics vary from study to study, read in as many as possible
  omics_sets <- strsplit(data_sets, "_")
  # Omit clinical information or custom mapping character string split from the 'omics portion
  # JC : I am still not a fan of importing an already made map and would like to return this line back 
  omics_sets <- omics_sets[-which(unlist(lapply(omics_sets, FUN=function(z) { any(c("clinical", "map") %in% z) })))]
  # 'omics names is concatenated from all strsplit prior to last element (which is presumably "_STUDYNAME.{rda,RData}")
  omics_names <- unlist(lapply(omics_sets, FUN=function(z) { paste(z[-length(z)], collapse="_") }))
  # Reconstruct paths to omics sets
  omics_sets <- paste0("data-raw/", unlist(lapply(omics_sets, FUN=function(z) { paste(z, collapse="_") })))
  # Construct a list of 'get' objects for the various omics, load them into environment then construct a list
  omics <- list()
  for(f in omics_sets){
    omic <- load(f)
    omics[[length(omics) + 1]] <- get(omic)
  }
  # Give correct name for each omics list member
  names(omics) <- omics_names
  if(verb == TRUE) print(paste0("Omics: ", paste(omics_names, collapse = ", ")))
  rm(f, omic, omics_names, omics_sets)
  
  ## Construct whole path to data sets
  data_sets <- paste0("data-raw/", data_sets)  
  
  pheno_name <- load(grep("data-raw/clinical_.*.RData", data_sets, value=TRUE))
  pheno_object <- get(pheno_name) 
  
  map <- data.frame(assay = character(0), primary = character(0), colname = character(0))
  for(i in 1:length(omics)){
    map <- rbind(map, data.frame(assay = names(omics)[[i]], 
                                 colname = colnames(omics[[i]]))
    )
  }
  
  sample_num <- stringr::str_count(pheno_object$sample_name[[1]], "\\|") + 1
  sample_num_names <- paste0("X",1:sample_num)
  
  map <- map %>% 
    dplyr::left_join(pheno_object %>% 
                       dplyr::select(primary = .data$patient_id, .data$sample_name) %>%
                       tidyr::separate(.data$sample_name, sample_num_names, 
                                       sep = "\\|") %>% 
                       dplyr::mutate_at(sample_num_names, 
                                        ~ gsub(".*: ", "", .)) %>% 
                       dplyr::mutate_at(sample_num_names, 
                                        ~dplyr::na_if(., "NA")) %>%
                       tidyr::pivot_longer(!.data$primary, names_to = "omic",
                                           values_to = "sample_name") %>% 
                       dplyr::select(-omic), 
                     by = c("colname" = "sample_name")) 
  
  # suggest we also pull the sample_name - MAE object doesnt want it here anyways 
  clinical_object <- pheno_object %>%
    dplyr::distinct(.data$patient_id, .keep_all = TRUE)
  
  clinical_object <- clinical_object %>% 
    dplyr::filter(.data$patient_id %in% map$primary)
  
  row.names(clinical_object) <- clinical_object$patient_id
  
  # Generate a MAE-object, generalization to various data compositions done above
  if(verb == TRUE) print("Final reformatting")  
  mae_object <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = MultiAssayExperiment::ExperimentList(omics),
    colData = as.data.frame(clinical_object),
    sampleMap = as.data.frame(map)
  )

  if(verb == TRUE) print(paste("MAE-object successfully created for", study_name))
  
  return(mae_object)
}

