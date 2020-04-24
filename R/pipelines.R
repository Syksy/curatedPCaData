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
	keys <- AnnotationDbi::mappedkeys(hgu133a.db::hgu133aGENENAME)
	nam <- names(as.character(hgu133a.db::hgu133aALIAS2PROBE)[match(rownames(GEX_Sun), as.character(hgu133a.db::hgu133aALIAS2PROBE))])
	nam[is.na(nam)] <- "NA"
	rownames(GEX_Sun) <- make.unique(nam)
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
	# Now returning ExpressionSet
	GEX_Sun	
}
