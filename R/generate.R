# Download gene expression from GEO using study specific id and process it
generate_gex_geo <- function(
  # Base directory for processing data, note a substantial amount of free HD space is required for raw data
  file_directory, 
  ## Allowed GEO ids:
  # "GSE21032" : Taylor et al.
  # "GSE25136" : Sun et al.
  geo_code = "GSE21032", # By default Taylor et al.
  # Whether downloaded and intermediate files ought to be cleaned up (deleted)
  cleanup = TRUE, 
  # Function to collapse probe(s) or select a probe, e.g. mean, median, or function that picks a probe with high variance
  collapseFUN = function(z) {apply(z, MARGIN = 2, FUN = stats::median)}, 
  # Function for cleaning rows/cols where GEO samples returned NaN or similar non-finite values only
  cleanFUN = janitor::remove_empty,
  ## Old parameter in Jordan's branch:
  #collapse_probes = function(z) {apply(z, MARGIN = 2, FUN = stats::median)}, # Function to collapse probe(s) or select a probe, e.g. mean, median, or function that picks a probe with high variance
  ...
){
  if(!missing(file_directory)) here::set_here(file_directory)
  # Supplementary files include the raw CEL files
  supfiles <- GEOquery::getGEOSuppFiles(geo_code)

  # Open the tarball(s)
  utils::untar(tarfile = rownames(supfiles))
  
  ##
  # Sun et al.
  ##
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
	gex <- do.call("rbind", by(as.matrix(affy::exprs(gex)), INDICES=nam, FUN=collapseFUN))
    
    ## TODO: Below snippet from Jordan, to be discussed; annotates new names
    #
    #compare_names <- data.frame(original = row.names(gex),
    #                            current = limma::alias2SymbolTable(row.names(gex),
    #                                                               species="Hs"))
    #duplicated_hugo_symbols <- compare_names[duplicated(compare_names$current),]$current
    #
    #compare_names <- compare_names %>% 
    #  dplyr::mutate(new_names = dplyr::case_when(
    #    current %in% duplicated_hugo_symbols ~ original,
    #    TRUE ~ current
    #  ))
    #
    #row.names(gex) <- compare_names$new_names
    #
    
	# Sort genes to alphabetic order for consistency
	gex <- gex[order(rownames(gex)),]
  }
  ##
  # Taylor et al.
  ##
  else if (geo_code == "GSE21032") {

	## Fetch supplementary files from Taylor et al. from GEO
	## Be wary as the tarball from Taylor et al. is 24.2 Gb 
	#supfiles <- GEOquery::getGEOSuppFiles('GSE21032')
	# Step in the internal directory created by GEOquery
	## TODO: Required in current pipeline?
	#setwd("GSE21032")
	# Open the tarball
	#utils::untar(tarfile=rownames(supfiles))
	
	# Read in the CEL files - note: requires a substantial amount of RAM for all 370 samples
	CELs <- oligo::read.celfiles(affy::list.celfiles())	
	
	# Perform RMA normalization
	RMAs <- oligo::rma(CELs)
	
	# Obtain gene and sample information
	featureData(RMAs) <- oligo::getNetAffx(RMAs, "transcript")
	# GSM######-type names from GEO
	nam0 <- unlist(lapply(strsplit(affy::list.celfiles(), "_"), FUN=function(z) z[[1]])) 
	# Two naming conventions if the files; picking the PCA###-style 
	nam1 <- unlist(lapply(strsplit(affy::list.celfiles(), "_"), FUN=function(z) z[[3]])) 
	nam2 <- gsub(".CEL.gz", "", unlist(lapply(strsplit(affy::list.celfiles(), "_"), FUN=function(z) z[[4]])))
	# Some samples were suffixed with HuEx, while others had Exonl prefix
	nam <- paste(nam0, "_", ifelse(nam1 == "Exon1", nam2, nam1), sep="")

	# Extract gene names
	genenames <- unlist(lapply(fData(RMAs)[,"geneassignment"], FUN=function(z) { strsplit(z, " // ")[[1]][2] }))

	# Transform into a matrix and remove empty gene names
	gex <- as.matrix(Biobase::exprs(RMAs))
	gex <- gex[-which(is.na(genenames)),]
	rownames(gex) <- genenames[-which(is.na(genenames))]
	# Map the coventionally used Taylor sample names instead of GEO codes 
	# Compatible with e.g. cBioPortal sample names
	# Give unique names with GSM#####.{PCA,PAN}##### combination for uniqueness
	colnames(gex) <- nam
	
	# Sort genes to alphabetic order for consistency
	gex <- gex[order(rownames(gex)),]
  }
  ##
  # Unknown GEO id, throw an R error
  ##
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
  # Return numeric matrix
  gex <- as.matrix(gex)
  gex <- gex %>% janitor::remove_empty(which = c("rows", "cols"))
  
  # clean names to match conventions? 
  # gex <- gex %>% janitor::clean_names()
  
}


# Download copy number variant data from GEO using study specific id and process it
generate_cna_geo <- function(
  # Base directory for processing data, note a substantial amount of free HD space is required for raw data
  file_directory, 
  ## Allowed GEO ids:
  # "GSE21035" : Taylor et al.
  # "GSE54691" : Hieronymus et al.
  geo_code = "GSE21035", # By default Taylor et al.
  # Whether downloaded and intermediate files ought to be cleaned up (deleted)
  cleanup = TRUE, 
  # Function for cleaning rows/cols where GEO samples returned NaN or similar non-finite values only
  cleanFUN = janitor::remove_empty,
  ...
){
  if(!missing(file_directory)) here::set_here(file_directory)
  # Supplementary files include the raw CEL files
  supfiles <- GEOquery::getGEOSuppFiles(geo_code)

  # Open the tarball(s)
  utils::untar(tarfile = rownames(supfiles))

  # Handle various GEO ids
  
  ##
  # same rCGH pipeline applied to datasets:
  # Taylor et al.
  # Hieronymus et al.
  ##
  if(geo_code %in% c("GSE21035", "GSE54691")){
  	# For now, the package 'rCGH' has to be available in the workspace, otherwise below functions will fail on e.g. rCGH::adjustSignal and when trying to find 'hg18'
  	library("rCGH")
  	# Read in Agilent 2-color data
	cna <- lapply(base::list.files(), FUN=function(z) { 
		try({
			cat("\n\nProcessing: ",z,"\n\n") 
			# Taylor et al.
			if(geo_code == "GSE21035"){
				rCGH::readAgilent(z, genome="hg38", sampleName=gsub(".txt.gz", "", z)) 
			# Hieronymus et al.
			}else if(geo_code == "GSE54691"){
				rCGH::readAgilent(z, genome="hg19", sampleName=gsub(".txt.gz", "", z)) 
			}
		})
	})
	#> list.files()[which(unlist(lapply(cna, FUN=class))=="try-error")]
	#[1] "GSM525755.txt" "GSM525763.txt"
	# Some files appear broken in Taylor et al; missing columns?
	
	# Omit data that could not be succcessfully read
	cna <- cna[-which(lapply(cna, FUN=class)=="try-error")]
	
	# Signal adjustments
	cna <- lapply(cna, FUN=function(z){
		try({
			rCGH::adjustSignal(z) 
		})
	})
	# Segmentation
	cna <- lapply(cna, FUN=function(z){
		try({
			rCGH::segmentCGH(z) 
		})
	})
	# EM-algorithm normalization
	cna <- lapply(cna, FUN=function(z){
		try({
			rCGH::EMnormalize(z) 
		})
	})
	# Remove additional suffixes from sample names
	cna <- lapply(cna, FUN=function(z){ 
		try({
			if(!"rCGH-Agilent" %in% class(z)){
				stop("This function is intended for Agilent aCGH analyzed with rCGH R Package (class \'rCGH-Agilent\')")		
			}
			# e.g. transform "GSM525575.txt|.gz" -> "GSM525575"
			z@info["sampleName"] <- gsub(pattern=".gz|.txt", replacement="", z@info["fileName"])
			z
		}) 
	})
	# Save sample names separately (of 'length(cna)')
	samplenames <- unlist(lapply(cna, FUN=function(z) { z@info["sampleName"] }))
	# Get segmentation table
	cna <- lapply(cna, FUN=function(z){
		try({
			rCGH::getSegTable(z)
		})
	})
	# Get per-gene table
	cna <- lapply(cna, FUN=function(z){
		try({
			rCGH::byGeneTable(z)
		})
	})
	# Extract all unique gene symbols present over all samples
	genenames <- unique(unlist(lapply(cna, FUN=function(z) { z$symbol })))
	# Bind genes to rows, name samples afterwards
	cna <- do.call("cbind", lapply(cna, FUN=function(z){
		# Return CNAs as Log2Ratios
		z[match(genenames, z$symbol), "Log2Ratio"]
	}))
	# Name rows and columns to genes and sample names, respectively
	rownames(cna) <- genenames
	colnames(cna) <- samplenames
	# CNA matrix is ready
	cna <- as.matrix(cna)
  }
  ##
  # Other (placeholder)
  ##
  else if(geo_code == ""){
  
  }
  ##
  # Unknown
  ##
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
  # Return numeric matrix
  as.matrix(cleanFUN(cna))
}


# Download generic 'omics data from cBioPortal using dataset specific query
# To examine available profiles after establishing mycgds-connection:
# cgdsr::getCaseLists(mycgds, cancerStudy="prad_tcga_pub")
# cgdsr::getGeneticProfiles(mycgds, cancerStudy="prad_tcga_pub")
generate_cbioportal <- function(
  # By default used the gene symbol from package data
  genes = sort(unique(curatedPCaData:::curatedPCaData_genes$hgnc_symbol)), # List of gene symbols to iterate over
  # geneticProfiles; Allowed platform/data combinations:
  # "prad_tcga_pub_rna_seq_v2_mrna" : TCGA GEX
  # "prad_tcga_pub_gistic" : TCGA CNA (GISTIC)
  # "prad_tcga_pub_linear_CNA" : TCGA CNA (Capped relative linear copy-number values)
  # "prad_mskcc_mrna_median_Zscores" : Taylor et al. GEX (z-score normalized)
  # "prad_mskcc_cna" : Taylor et al. CNA
  geneticProfiles = "prad_tcga_pub_rna_seq_v2_mrna", # for cgdsr calls, platform and dataset specific string
  # caseList; Allowed platform/data combinations:
  # "prad_tcga_pub_sequenced" : TCGA
  # "prad_mskcc_sequenced" : Taylor et al. 
  caseList = "prad_tcga_pub_sequenced", # for cgdsr calls, platform and dataset specific string
  # For Sys.sleep to not fetch too fast from cBio API
  delay = 0.05, 
  # Amount of genes fetched at a single API call - max 1000
  splitsize = 100, 
  # Function for cleaning rows/cols where cBio returned NaN or similar non-finite values only
  # clean_columns = janitor::remove_empty,
  # Verbosity
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
  genesplit <- rep(1:ceiling(length(genes)/splitsize), each=splitsize)[1:length(genes)]
  splitgenes <- split(genes, f=genesplit)
  # Fetch split gene name lists as separate calls
  pb <- progress::progress_bar$new(total = length(splitgenes))
  # Bind the API calls as per columns
  gex <- as.matrix(do.call("cbind", lapply(1:length(splitgenes), FUN=function(z){
      if(verb == TRUE) pb$tick()
      # Sleep if necessary to avoid API call overflow
      Sys.sleep(delay)
      # Fetch each split gene name list from the URL-based API, essentially a wrapper for cgdsr's own function
      cgdsr::getProfileData(mycgds, genes=splitgenes[[z]], geneticProfiles=geneticProfiles, caseList=caseList)
    })))  
  
  gex <- t(gex)
  
  gex <- gex %>% janitor::remove_empty(which = c("rows", "cols"))
  
}
