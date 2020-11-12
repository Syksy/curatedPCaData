#' Download gene expression from GEO using study specific id and process it
#' @param file_directory character string indicating path for downloading raw 
#' GEO data
#' @param geo_code character string indicating name of GEO dataset
#' @param cleanup logical value to remove intermediate files 
#' @param collapse_fun function to collapse probe(s) or select a probe, 
#' e.g. mean, median, or function that picks a probe with high variance
#' @param ... additional arguments
#' 
generate_gex_geo <- function(
  file_directory, 
  geo_code = c("GSE21032", # Taylor et al. TODO: Alternative more specific accession code "GSE21034" for GEX
               "GSE25136" # Sun et al.
               ), 
  cleanup = TRUE, 
  collapse_fun = function(z) {apply(z, MARGIN = 2, FUN = stats::median)},
  ...
){
  if(!missing(file_directory)) here::set_here(file_directory)
  # Supplementary files include the raw CEL files
  supfiles <- GEOquery::getGEOSuppFiles(geo_code)

  # Open the tarball(s)
  utils::untar(tarfile = rownames(supfiles))
  
  # Sun et al. -----
  
  if(geo_code == "GSE25136"){
	# Make sure to function in a working directory where the are no other tarballs present
	gz_files <- list.files()
	gz_files <- gz_files[grep(".gz", gz_files)]

	# Read Affymetrix MA
	Sun <- affy::ReadAffy()
	colnames(affy::exprs(Sun)) <- gsub(".gz|.CEL", "", colnames(Sun))

	# Careful not to mask 'rma' from 'affy' by the 'rma' from 'oligo'
	gex <- affy::rma(Sun)

	# Extracting .CEL and packaging names from the GEO-compatible sample names
	colnames(gex) <- gsub(".CEL.gz", "", colnames(affy::exprs(gex)))

	# Find gene annotations
	keys <- AnnotationDbi::mappedkeys(hgu133a.db::hgu133aGENENAME)
	nam <- names(as.character(hgu133a.db::hgu133aALIAS2PROBE)[match(rownames(gex),
								    as.character(hgu133a.db::hgu133aALIAS2PROBE))])
	nam[is.na(nam)] <- "NA"
	# Collapse probes
	gex <- do.call("rbind", by(as.matrix(affy::exprs(gex)), INDICES = nam, 
	                           FUN = collapse_fun))
	
    # Sun et al does not use Hugo names so the below code makes the names for sun 
	  # et al comparable with those in tcga/taylor
    compare_names <- data.frame(original = row.names(gex),
                               current = limma::alias2SymbolTable(row.names(gex),
                                                                  species="Hs"))
    duplicated_hugo_symbols <- compare_names[duplicated(compare_names$current),]$current

    compare_names <- compare_names %>%
     dplyr::mutate(new_names = dplyr::case_when(
       current %in% duplicated_hugo_symbols ~ original,
       TRUE ~ current
     ))

    row.names(gex) <- compare_names$new_names
    
	# Sort genes to alphabetic order for consistency
	gex <- gex[order(rownames(gex)),]
  }

  # Taylor et al.-----
  
  else if (geo_code == "GSE21032") { # TODO: Alternative more specific accession code "GSE21034"
	
	# Read in the CEL files - note: requires a substantial amount of RAM for all 370 samples
	CELs <- oligo::read.celfiles(affy::list.celfiles())	
	
	# Perform RMA normalization
	RMAs <- oligo::rma(CELs)
	
	# Obtain gene and sample information
	Biobase::featureData(RMAs) <- oligo::getNetAffx(RMAs, "transcript")
	# GSM######-type names from GEO
	nam0 <- unlist(lapply(strsplit(affy::list.celfiles(), "_"), 
	                      FUN = function(z) z[[1]])) 
	# Two naming conventions for the files; picking the PCA###-style 
	nam1 <- unlist(lapply(strsplit(affy::list.celfiles(), "_"), 
	                      FUN = function(z) z[[3]])) 
	nam2 <- gsub(".CEL.gz", "", unlist(lapply(strsplit(affy::list.celfiles(), "_"),
	                                          FUN = function(z) z[[4]])))
	# Some samples were suffixed with HuEx, while others had Exonl prefix
	nam <- paste(nam0, "_", ifelse(nam1 == "Exon1", nam2, nam1), sep="")

	# Extract gene names
	genenames <- unlist(lapply(Biobase::fData(RMAs)[,"geneassignment"], 
	                           FUN = function(z) { strsplit(z, " // ")[[1]][2] }))

	# Transform into a matrix and remove empty gene names
	gex <- as.matrix(Biobase::exprs(RMAs))
	gex <- gex[-which(is.na(genenames)),]
	rownames(gex) <- genenames[-which(is.na(genenames))]
	# Map the coventionally used Taylor sample names instead of GEO codes 
	# Compatible with e.g. cBioPortal sample names
	# Give unique names with GSM#####.{PCA,PAN}##### combination for uniqueness
	##colnames(gex) <- nam
	# Use the unique GSM###-names
	colnames(gex) <- unlist(lapply(nam, FUN=function(z){ strsplit(z, "_")[[1]][1] }))
	
	# Sort genes to alphabetic order for consistency
	gex <- gex[order(rownames(gex)),]
  }

  # Unknown GEO id (throw an error) -----

  else{
  	stop("Unknown GEO id, see allowed parameter values for geo_code")
  }

  # Remove downloaded files
  if(cleanup){
    # First GEO download
    file.remove(rownames(supfiles))
    # Tarballs
    file.remove(gz_files)
    # Remove empty folder
    file.remove(paste0(here::here(), "/", geo_code))
  }

  gex <- as.matrix(gex)
  gex <- gex %>% janitor::remove_empty(which = c("rows", "cols"))
  
}


#' Download copy number variant data from GEO using study specific id and process it
#' 
#' @param file_directory character string indicating path for downloading raw 
#' GEO data
#' @param geo_code character string indicating name of GEO dataset. Default is "GSE21035"
#' (Taylor et al)
#' @param cleanup logical value to remove intermediate files 
#' @param ... additional arguments
#' 
generate_cna_geo <- function(
  file_directory, 
  geo_code = c("GSE21035", # Taylor et al.
               "GSE54691" # Hieronymus et al.
               ),
  cleanup = TRUE, 
  ...
){
  if(!missing(file_directory)) here::set_here(file_directory)
  # Supplementary files include the raw CEL files
  supfiles <- GEOquery::getGEOSuppFiles(geo_code)

  # Hieronymus et al. requires going to a subfolder
  # you should NEVER use setwd() in a package - causes massive downstream problems 
  # need to remove 
  # if(geo_code == "GSE54691"){
  #   setwd("GSE54691") # Contains both CNA input as well as clinical data txt.gz
  # }
  
  # Open the tarball(s)
  utils::untar(tarfile = rownames(supfiles))
  
  ##
  # same rCGH pipeline applied to datasets:
  # Taylor et al.
  # Hieronymus et al.
  ##
  if(geo_code %in% c("GSE21035", "GSE54691")){
  	# For now, the package 'rCGH' has to be available in the workspace,
    # otherwise below functions will fail on e.g. rCGH::adjustSignal and when trying to find 'hg18'
  	# Read in Agilent 2-color data
  	cna <- lapply(list.files(geo_code), FUN = function(z) { 
  		try({
  			cat("\n\nProcessing: ",z,"\n\n") 
  			# Taylor et al. (and & Hieronymus et al.)
  			#if(geo_code %in% c("GSE21035", "GSE54691")){
			rCGH::readAgilent(z, genome = "hg38", sampleName = gsub(".txt.gz", "", z)) 
  			# Hieronymus et al. - any reason to use hg19 since the locations are known?
  			#}else if(geo_code == "GSE54691"){
  			#	rCGH::readAgilent(z, genome = "hg19", sampleName = gsub(".txt.gz", "", z)) 
  			#}
  		})
  	})
	#> list.files()[which(unlist(lapply(cna, FUN=class))=="try-error")]
	#[1] "GSM525755.txt" "GSM525763.txt"
	# Some files appear broken in Taylor et al; missing columns?
	
	# Not all files can always be successfully processed
	tryerr <- which(lapply(cna, FUN = class) == "try-error")
	if(length(tryerr)>0){
		# Warn of try-errors
		warning(paste("Error while processing files: ", 
			# Collapse file names
			paste(list.files()[which(unlist(lapply(cna, FUN=class))=="try-error")], collapse=", ")
		))
		# Omit data that could not be succcessfully read
		cna <- cna[-tryerr]
	}
	
	# Signal adjustments
	cna <- lapply(cna, FUN = function(z){
		try({
			rCGH::adjustSignal(z) 
		})
	})
	# Segmentation
	cna <- lapply(cna, FUN = function(z){
		try({
			rCGH::segmentCGH(z) 
		})
	})
	# EM-algorithm normalization
	cna <- lapply(cna, FUN = function(z){
		try({
			rCGH::EMnormalize(z) 
		})
	})
	# Remove additional suffixes from sample names
	cna <- lapply(cna, FUN = function(z){ 
		try({
			if(!"rCGH-Agilent" %in% class(z)){
				stop("This function is intended for Agilent aCGH analyzed with rCGH R Package (class \'rCGH-Agilent\')")		
			}
			# e.g. transform "GSM525575.txt|.gz" -> "GSM525575"
			z@info["sampleName"] <- gsub(pattern = ".gz|.txt", replacement = "", z@info["fileName"])
			z
		}) 
	})
	# Save sample names separately (of 'length(cna)')
	samplenames <- unlist(lapply(cna, FUN = function(z) { z@info["sampleName"] }))
	# Reformat samplenames in Hieronymus
	if(geo_code == "GSE54691"){
		# Pick GSM-part in GSM###_PCA###
		samplenames <- unlist(lapply(samplenames, FUN=function(z) { strsplit(z, "_")[[1]][[1]]}))
	}	
	# Get segmentation table
	cna <- lapply(cna, FUN = function(z){
		try({
			rCGH::getSegTable(z)
		})
	})
	# Get per-gene table
	cna <- lapply(cna, FUN = function(z){
		try({
			rCGH::byGeneTable(z)
		})
	})
	# Extract all unique gene symbols present over all samples
	genenames <- unique(unlist(lapply(cna, FUN = function(z) { z$symbol })))
	# Bind genes to rows, name samples afterwards
	cna <- do.call("cbind", lapply(cna, FUN = function(z){
		# Return CNAs as Log2Ratios
		z[match(genenames, z$symbol), "Log2Ratio"]
	}))
	# Name rows and columns to genes and sample names, respectively
	rownames(cna) <- genenames
	colnames(cna) <- samplenames
	# CNA matrix is ready
	cna <- as.matrix(cna)
  }
  # Other (placeholder) - should probably place this in the begining so it 
  # doesnt try to run through all the beginning steps 
  else if(geo_code == ""){
    stop("Must supply a GEO id")
  }
  # Unknown
  else{
    stop("Unknown GEO id, see allowed parameter values for geo_code")
  }

  # Remove downloaded files
  # TODO: Make sure appropriate permissions exist for removing files
  if(cleanup){
    # First GEO download
    file.remove(rownames(supfiles))
    # Tarballs
    # TODO: Neither CNA pipeline produces .gz-files
    #file.remove(gz_files)
    # Remove empty folder
    file.remove(paste0(here::here(), "/", geo_code))
  }
  # Return numeric matrix
  cna <- as.matrix(cna)
  cna <- cna %>% janitor::remove_empty(which = c("rows", "cols"))
}


#' Download generic 'omics data from cBioPortal using dataset specific query
#' 
#' @param genes character vector of genes to query
#' @param geneticProfiles charatcer string of cbioportal genetic profiles
#' @param caseList charcter string of patient IDs for that genetic profile
#' @param delay numberic value for delay time between querying gene sets
#' @param splitsize number of genes in each query
#' @param verb logical value for displaying progress bar 
#' @param ... any additional arguments
#' 
generate_cbioportal <- function(
  genes = sort(unique(curatedPCaData_genes$hgnc_symbol)), 
  geneticProfiles = c("prad_tcga_pub_rna_seq_v2_mrna", #TCGA GEX 
                      "prad_tcga_pub_gistic", # TCGA CNA (GISTIC)
                      "prad_tcga_pub_linear_CNA", # TCGA CNA (Capped relative linear copy-number values)
                      "prad_mskcc_mrna_median_Zscores", # Taylor et al. GEX (z-score normalized)
                      "prad_mskcc_cna" # Taylor et al. CNA
                      ), # for cgdsr calls, platform and dataset specific string
  caseList = c("prad_tcga_pub_sequenced", # TCGA
               "prad_mskcc_sequenced" # Taylor et al. 
               ), # for cgdsr calls, platform and dataset specific string
  delay = 0.05, 
  splitsize = 100, 
  verb = TRUE,
  ...
){
  # If given genes is a list (with slots for various annotation types), try to extract hugo gene symbols
  if(class(genes)=="list"){
    genes <- genes$hgnc_symbol
  }
  # Establisigh connection to cBioPortal
  mycgds <- cgdsr::CGDS("http://www.cbioportal.org/")
  # Split gene name vector into suitable lengths
  genesplit <- rep(1:ceiling(length(genes)/splitsize), 
                   each = splitsize)[1:length(genes)]
  splitgenes <- split(genes, f = genesplit)
  # Fetch split gene name lists as separate calls
  pb <- progress::progress_bar$new(total = length(splitgenes))
  # Bind the API calls as per columns
  gex <- as.matrix(do.call("cbind", lapply(1:length(splitgenes), FUN = function(z){
      if(verb == TRUE) pb$tick()
      # Sleep if necessary to avoid API call overflow
      Sys.sleep(delay)
      # Fetch each split gene name list from the URL-based API, essentially a wrapper for cgdsr's own function
      cgdsr::getProfileData(mycgds, genes = splitgenes[[z]], 
                            geneticProfiles = geneticProfiles, caseList = caseList)
    })))  
  
  gex <- t(gex)
  
  gex <- gex %>% janitor::remove_empty(which = c("rows", "cols"))
  
}
