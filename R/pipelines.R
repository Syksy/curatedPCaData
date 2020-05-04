#' Generate gene expression data for Sun et al. 
#' 
#' @param file_directory Character string indicating directory for downloading files. 
#' Files may be large, so please check that there is sufficient free space. 
#' If NULL then files are downloaded into current directory.
#' @param cleanup Logical. Remove tarballs and other files from working directory.
#' @param ... Additional arguments. 
#' @return Gene expression object of a particular type 
#' @example 
#' GEX_Sun <- Generate_GEX_Sun()
Generate_GEX_Sun <- function(
	file_directory, 
	cleanup = FALSE, 
	collapseFUN = function(z) { apply(z, MARGIN=2, FUN=median) }, # Function to collapse probe(s) or select a probe, e.g. mean, median, or function that picks a probe with high variance
	...
){
	if(!missing(file_directory)) setwd(file_directory) # exchange setwd with here::here()
	# Supplementary files include the raw CEL files
	supfiles <- GEOquery::getGEOSuppFiles('GSE25136')
	# Download size: 269.5 MB
	# Open the tarball(s)
	utils::untar(tarfile=rownames(supfiles))
	# Make sure to function in a working directory where the are no other tarballs present
	supfiles2 <- list.files()
	supfiles2 <- supfiles2[grep(".gz", supfiles2)]
	#invisible(lapply(supfiles2, FUN=gunzip))
	# Read Affymetrix MA
	Sun <- affy::ReadAffy()
	colnames(affy::exprs(Sun)) <- gsub(".gz|.CEL", "", colnames(Sun))
	# Careful not to mask 'rma' from 'affy' by the 'rma' from 'oligo'
	GEX_Sun <- affy::rma(Sun)
	# Removing .CEL and packaging names from the GEO-compatible sample names
	colnames(GEX_Sun) <- gsub(".CEL.gz", "", colnames(affy::exprs(GEX_Sun)))
	#GEX_Sun <- Sun2
	# keys <- AnnotationDbi::mappedkeys(hgu133a.db::hgu133aGENENAME)
	nam <- names(as.character(hgu133a.db::hgu133aALIAS2PROBE)[match(rownames(GEX_Sun), as.character(hgu133a.db::hgu133aALIAS2PROBE))])
	nam[is.na(nam)] <- "NA"
	GEX_Sun <- do.call("rbind", by(as.matrix(exprs(GEX_Sun)), INDICES=nam, FUN=collapseFUN))
	#rownames(GEX_Sun) <- make.unique(nam)
	# Remove downloaded files
	if(cleanup){
		# First GEO download
		file.remove(rownames(supfiles))
		# Tarballs
		file.remove(supfiles2)
		# Remove empty folder
		file.remove(paste0(getwd(),"/GSE25136"))
	}
	# TODO: Transform into a MultiAssayExperiment-object prior to returning object (MAE_Sun)
	# Return numeric matrix
	as.matrix(GEX_Sun)
}
<<<<<<< HEAD
=======

Generate_GEX_TCGA <- function(
	genes, # List of gene symbols to iterate over
	...
){
	#http://www.cbioportal.org/study?id=prad_tcga#summary
	#mycgds <- cgdsr::CGDS("http://www.cbioportal.org/public-portal/")
	mycgds <- cgdsr::CGDS("http://www.cbioportal.org/")
	# Get the summary lists for 'omics and case profiles	
	TCGA <- cgdsr::getCaseLists(mycgds,"prad_tcga")
	TCGA_genetic <- cgdsr::getGeneticProfiles(mycgds,"prad_tcga")
	# Use the wrapper function to iterative calls for split gene lists
	#GEX_TCGA <- curatedPCaData::getProfileDataWrapper(
	GEX_TCGA <- getProfileDataWrapper(
		x=mycgds, # cgdsr object
		#genes=curatedPCaData:::.getGeneNames()$hgnc, # All unique gene symbols
		genes=genes, # All unique gene symbols
		geneticProfiles="prad_tcga_rna_seq_v2_mrna", # mRNA expression
		caseList="prad_tcga_sequenced", # Case list
		verb = 2
	)
	# TODO: Filter out low quality samples based on list provided by Travis
	# Return GEX for TCGA
	GEX_TCGA
}

Generate_cBioPortal <- function(
	genes, # List of gene symbols to iterate over
	geneticProfiles, # for cgdsr calls, platform and dataset specific string
	caseList, # for cgdsr calls, platform and dataset specific string
	...
){
	#http://www.cbioportal.org/study?id=prad_tcga#summary
	mycgds <- cgdsr::CGDS("http://www.cbioportal.org/")
	# Use the wrapper function to iterative calls for split gene lists
	dat <- getProfileDataWrapper(
		x=mycgds, # cgdsr object
		#genes=curatedPCaData:::.getGeneNames()$hgnc, # All unique gene symbols
		genes=genes, # All unique gene symbols
		geneticProfiles="prad_tcga_rna_seq_v2_mrna", # Omics profile
		caseList="prad_tcga_sequenced" # Case list
	)
	# TODO: Filter out low quality samples based on list provided by Travis
	# Return omics matrix
	dat
}

####
#
# Supporting variables
#
####

# Extract whole list of desired gene names in various formats
genes <- .getGeneNames()

####
#
# TCGA 
# GEX + CNA
#
####

## Get the summary lists for 'omics and case profiles (gives names for omics 'geneticProfiles' and 'caseList')
# cgdsr::getCaseLists(mycgds,"prad_tcga")
# cgdsr::getGeneticProfiles(mycgds,"prad_tcga")

# Running TCGA GEX:
GEX_TCGA <- Generate_cBioPortal(genes = genes$hgnc, geneticProfiles="prad_tcga_rna_seq_v2_mrna", caseList="prad_tcga_sequenced")

# Running TCGA CNA:
CNA_TCGA <- Generate_cBioPortal(genes = genes$hgnc, geneticProfiles="prad_tcga_gistic", caseList="prad_tcga_sequenced")

####
#
# Sun, et al.
# GEX
#
####

# Running Sun, et al.:
GEX_Sun <- Generate_GEX_Sun()

>>>>>>> origin/master
