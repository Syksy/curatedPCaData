#' Download gene expression from GEO using study specific id
generate_gex_geo <- function(
  file_directory, 
  ## Allowed GEO ids:
  # "GSE25136" : Sun et al.
  # "GSE21032" : Taylor et al.
  geo_code = "GSE25136", # code for Sun et al. (Taylor et al. - GSE21032)
  cleanup = TRUE, 
  collapseFUN = function(z) {apply(z, MARGIN = 2, FUN = stats::median)}, # Function to collapse probe(s) or select a probe, e.g. mean, median, or function that picks a probe with high variance
  # Function for cleaning rows/cols where GEO samples returned NaN or similar non-finite values only
  cleanFUN = janitor::remove_empty,
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
	nam <- paste(nam0, ".", ifelse(nam1 == "Exon1", nam2, nam1), sep="")

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
  as.matrix(cleanFUN(gex))
}

#' Download generic 'omics data from cBioPortal using dataset specific query
generate_cbioportal <- function(
  # By default used the gene symbol from package data
  genes = sort(curatedPCaData:::tcga_gene_names$hgnc), # List of gene symbols to iterate over
  # geneticProfiles; Allowed platform/data combinations:
  # "prad_tcga_pub_rna_seq_v2_mrna" : TCGA GEX
  # "prad_tcga_pub_gistic" : TCGA CNA (GISTIC)
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
  cleanFUN = janitor::remove_empty,
  # Verbosity
  verb = TRUE,
  ...
){
  mycgds <- cgdsr::CGDS("http://www.cbioportal.org/")
  genesplit <- rep(1:ceiling(length(genes)/splitsize), each=splitsize)[1:length(genes)]
  #genesplit <- rep(1:ceiling(length(genes$hgnc)/splitsize), each=splitsize)[1:length(genes$hgnc)]
  splitgenes <- split(genes, f=genesplit)
  #splitgenes <- split(genes$hgnc, f=genesplit)
  # Fetch split gene name lists as separate calls
  pb <- progress::progress_bar$new(total = length(splitgenes))
  # Transpose to have genes as row and samples as columns
  cleanFUN(t(
    # Bind the API calls as per columns
    as.matrix(do.call("cbind", lapply(1:length(splitgenes), FUN=function(z){
      # if(verb == TRUE) print(paste("Downloading gene set", z, "of", length(splitgenes), "from cBioportal"))
      if(verb == TRUE) pb$tick()
      # Sleep if necessary to avoid API call overflow
      Sys.sleep(delay)
      # Fetch each split gene name list from the URL-based API, essentially a wrapper for cgdsr's own function
      cgdsr::getProfileData(mycgds, genes=splitgenes[[z]], geneticProfiles=geneticProfiles, caseList=caseList)
    }))
  ))  
}
