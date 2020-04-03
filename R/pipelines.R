

Generate_GEX_Sun <- function(
	wd, # Working directory; may be needed to download large files, so ought to be one with sufficient free space
	cleanup = FALSE, # Remove tarballs etc in working directory; careful not to remove any other potential *tz, trying to be selective, so turn on manually
	...
){
	if(!missing(wd)) setwd(wd)
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
	colnames(exprs(Sun)) <- gsub(".gz|.CEL", "", colnames(Sun))
	# Careful not to mask 'rma' from 'affy' by the 'rma' from 'oligo'
	GEX_Sun <- affy::rma(Sun)
	# Removing .CEL and packaging names from the GEO-compatible sample names
	colnames(GEX_Sun) <- gsub(".CEL.gz", "", colnames(exprs(GEX_Sun)))
	#GEX_Sun <- Sun2
	keys <- mappedkeys(hgu133a.db::hgu133aGENENAME)
	nam <- names(as.character(hgu133a.db::hgu133aALIAS2PROBE)[match(rownames(GEX_Sun), as.character(hgu133a.db::hgu133aALIAS2PROBE))])
	nam[is.na(nam)] <- "NA"
	rownames(GEX_Sun) <- make.unique(nam)
	# Remove downloaded files
	if(cleanup){
		# First GEO download
		file.remove(rownames(supfiles))
		# Tarballs
		file.remove(supfiles2)
	}
	# TODO: Transform into a MultiAssayExperiment-object prior to returning object (MAE_Sun)
	# Now returning ExpressionSet
	GEX_Sun	
}

# Running Sun, et al.:
GEX_Sun <- Generate_GEX_Sun()


