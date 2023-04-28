###
#
# getPCa:
# Main function and additional support functions for fetching Prostate Cancer (PCa) 
# datasets from ExperimentHub curated via curatedPCaData
#
###

.conditionToIndex <- function(startVec, testVec, FUN) {
  logmat <- vapply(startVec, FUN, logical(length(testVec)))
  apply(logmat, 1L, any)
}

.getResources <- function(ExperimentHub, resTable, verbose) {
  fileNames <- stats::setNames(resTable[["RDataPath"]], resTable[["Title"]])
  resources <- lapply(resTable[["Title"]], function(res) {
    if (verbose)
      message("Working on: ", gsub("\\.rds", "", basename(res)))
    query(ExperimentHub, res)[[1L]]
 
  })
  
  names(resources) <- resTable[["Title"]]
  resources
}

.test_eh <- function(...) {
  tryCatch({
    ExperimentHub(...)
  }, error = function(e) {
    emsg <- conditionMessage(e)
    if (grepl("Timeout", emsg))
      warning("[experimenthub.bioconductor.org] timeout, localHub=TRUE",
              call.=FALSE)
    ExperimentHub(..., localHub = TRUE)
  })
}

#' Create a MultiAssayExperiment from specific assays and cohorts
#'
#' @description curatedPCaData provides \linkS4class{MultiAssayExperiment} 
#' container objects that are constructed from ExperimentHub.
#' User provides PCa data set name (see list XXXX) and YYYY.
#'
#' @details This function will check against available resources in
#' ExperimentHub. 
#' For a list of datasets, see the XXXX
#'
#' @param dataset character() of PCa cancer cohort names
#'     (e.g., 'abida')
#'
#' @param slots character() A vector of PCa assays. If not included, returns all
#'     available for the selected dataset;
#'     see below for more details
#'
#' @param timestamp character(1) "20230215" indicating the data version to obtain 
#'     from `ExperimentHub`. See `version` section details. XXXXXXXXXXXXXXX
#'
#' @param ... Additional arguments passed on to the
#'     \code{\link[ExperimentHub:ExperimentHub-class]{ExperimentHub}}
#'     constructor
#'
#' @param verbose logical(1) Whether to show the dataset currenlty being
#'     (down)loaded (default TRUE)
#'
#' @section Available Assays:
#'
#' Thelist of ExperimentList assay names and their descriptions.
#' These assays can be entered as part of the \code{slots} argument in the 
#' main function. 
#' \preformatted{
#'
#' ExperimentList data types   Description
#' ----------------------------------------------------------------------------
#' 
#'   gex.rma            Gene expression values
#'   gex.logq              	   
#'   gex.relz			
#'   gex.logr
#'   gex.rsem.log
#'   cna.gistic         Copy number alteration
#'   cna.logr
#'   mut                Somatic mutations
#'   cibersort				
#'   xcell
#'   epic
#'   quantiseq
#'   mcp
#'   estimate
#'   scores      
#' }
#'
#' @section timestamp: XXXXXX
#' The timestamp is updated in case the data is updated. In this case, this 
#' section describes the changes made for new timestamps.
#'
#' @seealso curatedPCaData-package
#'
#' @return a \linkS4class{MultiAssayExperiment} of the specified assays and
#' cancer code
#'
#' @examples
#'
#' mae_taylor <- getPCa(
#'     dataset = "taylor", timestamp = "20230215"
#' )
#'
#'
#' @md
#'
#' @export getPCa
getPCa <- function(
	# Dataset name
	dataset, 
	# Data slots to retrieve (i.e. user can subset to just desired data)
	slots, 
	# Timestamps of data from ExperimentHub; allowed values: '20230215'
	timestamp,
	# Verbosity
	verbose = FALSE, 
	# Additional parameters
	...
)
{
	if(missing(dataset)){
		stop("Select dataset; see ?curatedPCaData")
	}

	if(length(dataset)>1){
		stop("Select only one dataset at a time.")
	}

	dataset <- tolower(dataset)  # datasets are saved as lower case

	## Change this if new update (with new timestamp) is added
	if(missing(timestamp)){
		timestamp <- "20230215"
	}

	# Update this if more timestamps becomes available
	if (any(!timestamp %in% c("20230215"))){
	  stop("'timestamp' contains a timestamp that is not '20230215'; see '?curatedPCaData'")}

	assays_file <- system.file("extdata", "metadata.csv", package = "curatedPCaData", mustWork = TRUE)
	assay_metadat <- read.csv(assays_file, stringsAsFactors = FALSE)

	# Separate names Title in assay_metadat to get dataset, assay and timestamp separately
	eh_assays <- assay_metadat[["Title"]]
	eh_assays_sep <- t(as.data.frame(strsplit(eh_assays,"_")))

	# Get only requested dataset
	dataId <- which(eh_assays_sep[,1]==dataset)

	if(length(dataId)==0){
	  stop("Dataset name is incorrect or not available.")
	}

	eh_assays_sep <- eh_assays_sep[dataId,]

	assaysAvail <- unique(eh_assays_sep[,2])  # Get available assays for selected dataset

	# Select user specified assays
	if(!missing(slots)){ # If nothing specific requested, return all
	   if(any(!slots %in% assaysAvail)){ # If user asks for something that is not available, 
		  stop(paste0(c("At least one of asked slots is not available. The available slots for this dataset are:",assaysAvail),collapse="  "))
	   }else{ # Select only requested assays
		  assaysAvail <- unique(c(slots,"colData","sampleMap"))
	   }
	}

	# Select assays by timestamp request. If more versions are added this has to be updated.
	# sampleMap is always selected by the latest timestamp given. If not available, latest is used.
	codeAssay <- c()
	latest <- as.character(max(as.numeric(timestamp)))
	
	for(i in assaysAvail){
		assayId <- which(eh_assays_sep[,2]==i)
		availableTimestamp <- eh_assays_sep[assayId,3]

		if(i == "sampleMap"){  # Use latest given timestamp for sampleMap, if not available, use latest in total
		  if(latest %in% availableTimestamp){
			selectedTimestamp <- latest
		  }else{
			selectedTimestamp <- as.character(max(as.numeric(availableTimestamp)))
		  }
		}else{
		  # Go through the user given timestamp request and get the first match, if no match, give latest
		  nomatch <- FALSE
		  n <- 1
		  while(!nomatch){
			# Select the first match, if no match, select latest
			if(timestamp[n] %in% availableTimestamp){
			  nomatch <- TRUE
			  selectedTimestamp <- timestamp[n]
			}else{
			  n<-n+1
			  if(n>length(timestamp)){
				nomatch <- TRUE
				selectedTimestamp <- as.character(max(as.numeric(availableTimestamp)))
			  }
			}
		  }
		}

		codeAssay <- c(codeAssay,paste0(dataset,"_",i,"_",selectedTimestamp))
	}
	
	# Get indices of selected dataset_assay_timestamp combos
	fileIdx <- .conditionToIndex(codeAssay, eh_assays, function(x) startsWith(eh_assays, x))
	fileMatches <- assay_metadat[fileIdx, c("Title", "DispatchClass")]

	# Sanitycheck
	if (!length(nrow(fileMatches))){
		stop("Cancer and data type combination(s) not available")
	}

	eh <- .test_eh()

	# Get the data
	assay_list <- .getResources(
		eh, assay_metadat[fileIdx, c("Title", "RDataPath")], verbose
	)

	# Get experiments and name them only by the assay/experiment/etc name
	# Omit colData and sampleMap since they are added separately
	cD_idx <- which(grepl("colData", names(assay_list), fixed = TRUE))
	sM_idx <- which(grepl("sampleMap", names(assay_list), fixed = TRUE))

	# Final sanity check
	if(length(cD_idx)!=1 | length(sM_idx)!=1){
		stop("Either colData or sampleMap is missing for this dataset and the MAE cannot be constructed, or there are multiple with same timestamp.")
	}

	eh_experiments <- ExperimentList(assay_list[-c(cD_idx,sM_idx)])
	names(eh_experiments) <- gsub("(^[a-z]*)_(.*)_(.*)", "\\2", names(eh_experiments))

	# Inform user
	cat(paste0("\nConstructing MultiAssayExperiment for study ", dataset, " from local cache downloaded from ExperimentHub"))

	# Return MAE
	MultiAssayExperiment::MultiAssayExperiment(
		experiments = eh_experiments,
		colData = assay_list[cD_idx][[1]],
		sampleMap = assay_list[sM_idx][[1]]
	)	
}

