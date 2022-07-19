#' Create multi-assay experiment object. 
#' 
#' @param study_name study identifier.
#' @param verb print output as MAE is being created 
#' @param ... additional arguments 
#' @return An MAE object containing clinical and multi-'omics data
#' 
#' @noRd
#' @keywords internal
create_mae <- function(
  # Valid study_name-parameters: <list latest valid study names>
  study_name,
  # Level of verbosity
  verb = TRUE,
  ...
){
  # Separate loading of RaggedExperiment required to be able to create MultiAssayExperiment::ExperimentList containing RaggedExperiment-objects
  requireNamespace("RaggedExperiment")
  # NOTE: Normally loading of packages should not be required
  # A better fix is probably possible, although extensive testing didn't provide one that would've helped 
  # with the 'Error values must be length 1 but FUN(X[[1]]) result is length 0 ...' error
  # when casting to MultiAssayExperiment::ExperimentList with the list containing both matrices and RaggedExperiment
  
  if(verb) print(paste("Starting to process:", study_name))
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
  if(verb) print(paste0("Omics: ", paste(omics_names, collapse = ", ")))
  #rm(f, omic, omics_names, omics_sets)
  
  ## Construct whole path to data sets
  data_sets <- paste0("data-raw/", data_sets)  
  
  pheno_name <- load(grep("data-raw/clinical_.*.RData", data_sets, value=TRUE))
  pheno_object <- get(pheno_name) 
  
  map <- data.frame(assay = character(0), primary = character(0), colname = character(0))
  for(i in 1:length(omics)){
  	# GEX/CNA stored as matrices; their sample names are column names
	if(any(class(omics[[i]]) %in% c("matrix","array"))){
	    map <- rbind(map, data.frame(assay = names(omics)[[i]], 
		colname = colnames(omics[[i]])
	      )
	    )
        # RaggedExperiments (for 'mut') store their sample names for mapping differently than ordinary matrices (for 'gex' and 'cna)
	}else{
	    map <- rbind(map, data.frame(assay = names(omics)[[i]], 
		colname = names(omics[[i]]@assays)
	      )
	    )
	}
  }
  
  sample_num <- stringr::str_count(pheno_object$sample_name[[1]], "\\|") + 1
  sample_num_names <- paste0("X",1:sample_num)
  
  # Format .data as NULL for visibility as it is an internal dplyr structure; will not pass R CMD check otherwise
  .data <- NULL
  
  map <- map |>
    dplyr::left_join(pheno_object |>
                       dplyr::select(primary = .data$patient_id, .data$sample_name) |>
                       tidyr::separate(.data$sample_name, sample_num_names, sep = "\\|") |>
                       dplyr::mutate_at(sample_num_names, ~ gsub(".*: ", "", .)) |>
                       dplyr::mutate_at(sample_num_names, ~ dplyr::na_if(., "NA")) |>
                       tidyr::pivot_longer(!.data$primary, names_to = "omic", values_to = "sample_name") |>
                       dplyr::select(-omic), 
                     by = c("colname" = "sample_name")) 
  
  clinical_object <- pheno_object |>
    dplyr::distinct(.data$patient_id, .keep_all = TRUE)
  
  clinical_object <- clinical_object |>
    dplyr::filter(.data$patient_id %in% map$primary)
  
  row.names(clinical_object) <- clinical_object$patient_id
  
  # Generate a MAE-object, generalization to various data compositions done above
  if(verb) print("Final reformatting")  
  mae_object <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = MultiAssayExperiment::ExperimentList(omics),
    colData = as.data.frame(clinical_object),
    sampleMap = as.data.frame(map)
  )

  if(verb) print(paste("Saving system date in attr 'mae_date' for sanity checking:", Sys.time()))
  attr(mae_object, 'mae_date') <- Sys.time()

  if(verb) print(paste("MAE-object successfully created for", study_name))
  
  return(mae_object)
}

