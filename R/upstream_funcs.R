####
##
## Upstream tools; hidden from export via .-notation for most part
##
####

###
#
# Fetching tools
#
###

#' Wrapper to fetch omics for larger gene lists from cBioPortal's URL-based interface and then bind them together
#'
#' @param x See parameter 'x' in cgdsr::getProfileData
#' @param genes See parameter 'genes' in cgdsr::getProfileData
#' @param geneticProfiles See parameter 'geneticProfiles' in cgdsr::getProfileData
#' @param caseList See parameter 'caseList' in cgdsr::getProfileData
#' @param delay Delay in seconds for Sys.sleep to not fetch too fast from cBioPortal API
#' @param splitsize How many genes / locations / IDs / etc are fetched at each call to cBioPortal API
#' @param verb Level of verbosity (integer); If call should be verbose; 0 = silent, 1 = info
#'
#' @return A matrix of corresponding omics-measurements, where the batch-wise fetches have been bound as columns.
#' 
#' @examples
#' # Fetching over all genes
#' # biomaRt-package for latest hg annotations etc
#' library(biomaRt)
#' ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
#' genes <- getBM(attributes=c('ensembl_gene_id', 'ensembl_transcript_id','hgnc_symbol','chromosome_name','start_position','end_position'), mart = ensembl)
#' # Using hgnc gene symbols 
#' genenames <- unique(genes$hgnc_symbol)
# # Using cBioPortal's package 'cgdsr' to access data from cBio
#' library(cgdsr)
#' # Create CGDS object
#' mycgds <- CGDS("http://www.cbioportal.org/public-portal/")
#' # http://www.cbioportal.org/study?id=prad_tcga#summary
#' # Get the summary lists for 'omics and case profiles
#' TCGA <- getCaseLists(mycgds,"prad_tcga")
#' TCGA_genetic <- getGeneticProfiles(mycgds,"prad_tcga")
#' # mRNA expression (RNA Seq V2 RSEM)
#' # NOTE: Below call(s) can take up a notable amount of time if full genomic profile is fetched.
#' 
#'\donttest{
#' TCGA_mRNA <- getProfileDataWrapper(x=mycgds, genes=genenames,geneticProfiles="prad_tcga_rna_seq_v2_mrna",caseList=
#'}
#'
#' @note Intended to work as a wrapper function for the function 'getProfileData' from package 'cgdsr'. Notice that most parameters are passed on to this underlying function as-is.
#' @author Teemu Daniel Laajala \email{teelaa@@utu.fi}
#' @export
getProfileDataWrapper <- function(
	x, 
	genes, 
	geneticProfiles, 
	caseList, 
	delay = 0.05, # For Sys.sleep to not fetch too fast from cBio API
	splitsize = 100, # How many genes are fetched at one time
	verb = 1 # If call should be verbose; 0 = silent, 1 = info
){
	genesplit <- rep(1:ceiling(length(genes)/splitsize), each=splitsize)[1:length(genes)]
	splitgenes <- split(genes, f=genesplit)
	# Fetch split gene name lists as separate calls
	as.matrix(do.call("cbind", lapply(1:length(splitgenes), FUN=function(z){
		if(verb>=1) print(paste("Calling", z, "of", length(splitgenes), "..."))
		# Sleep if necessary to avoid API call overflow
		Sys.sleep(delay)
		# Fetch each split gene name list from the URL-based API, essentially a wrapper for cgdsr's own function
		cgdsr::getProfileData(x=x, genes=splitgenes[[z]], geneticProfiles=geneticProfiles, caseList=caseList)
	})))	
}

#' Preliminary function for fetching gene names
.getGeneNames <- function(
	update = TRUE, # Whether Bioconductor packages require updating; sometimes a key package has been updated, but may also be optional
	ask = FALSE # Do not prompt user for 'a/s/n' package update selection by default
){
	# Bioconductor;
	# Warning! Old packages may cause some dependencies to fail, while their updating may fail if they're already in use for the session
	# Safest is to update all Bioconductor packages before analyses / running pipelines
	if (!requireNamespace("BiocManager", quietly = TRUE))
	    install.packages("BiocManager")
	BiocManager::install(version = "3.10", update=update, ask=ask)
	# BioConductor 'annotate' package for all sorts of conversions and genetic location info, etc
	# Annotation for microarrays
	# For Entrez <-> Hugo Gene Symbol mapping
	# http://bioconductor.org/packages/release/bioc/html/annotate.html
	if(!require(annotate)) {
		biocLite("annotate", suppressUpdates=T)
		library(annotate)
	}
	# Genome wide annotation for Human
	# For Entrez <-> Hugo Gene Symbol mapping (database)
	# https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html
	if(!require(org.Hs.eg.db)){
		biocLite("org.Hs.eg.db", suppressUpdates=T)
		library(org.Hs.eg.db)
	}

	# biomaRt
	if(!require(biomaRt)){
		biocLite("biomaRt", suppressUpdates=T)
		library(biomaRt)
	}

	# Fetch gene names
	ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
	genes <- getBM(attributes=c('ensembl_gene_id', 'ensembl_transcript_id','hgnc_symbol','chromosome_name','start_position','end_position'), mart = ensembl)
	# Using hgnc gene symbols 
	genenames <- unique(genes$hgnc_symbol)

	# Alternative approach using 'GenomicFeatures'
	if(!require(GenomicFeatures)){
		stop("Require 'GenomicFeatures' for GenomicRanges")
	}
	hg38_genes <- makeTxDbFromUCSC(genome="hg38", table="refGene")
	refseq_genes <- genes(hg38_genes)
	# Return 
	list(hgnc = genenames, hg38_genes = hg38_genes, refseq_genes = refseq_genes)
}

# Hidden function that renames missing sample names in rCGH files from the corresponding input file name
.renameSampleName <- function(
	x
){
	try({
		if(!"rCGH-Agilent" %in% class(x)){
			stop("This function is intended for Agilent aCGH analyzed with rCGH R Package (class \'rCGH-Agilent\')")		
		}
		# e.g. append "GSM525575.txt" -> ""GSM525575"
		x@info["sampleName"] <- gsub(pattern=".txt", replacement="", x=x@info["fileName"])
		x
	})
}

###
#
# Functions for the ICGC datasets
#
###

#' Function for downloading ICGC files from a certain release using http and then open the packaged tar.gz
.icgcDownload <- function(url){
	# Pick filename from the end of the URL
	filename <- strsplit(url, "/")
	filename <- filename[[1]][[length(filename[[1]])]]
	# Download file into parsed *.tsv.gz 
	download.file(url=url, destfile=filename)
	# gunzip the files open
	R.utils::gunzip(filename)
}

# Process exp_array.*.tsv and exp_seq.*.tsv from ICGC
.processICGCexp <- function(
	file,
	type = "array", # Should be either 'array' or 'seq' indicating if the input is microarray or RNA-seq
	verb = 0 # Level of verbosity as in every verb:th row print the location in the matrix that is being populated, 0 indicates no reporting of any kind by default
){
	# Read the file from a possible entry in current wd or have file location/name given as a parameter
	# Error handling done using 'stop' to prevent strange input files as this only works for ICGC raw downloads
	if(missing(file)){
		f <- list.files()
		if(type=="array"){
			file <- f[grep("exp_array", f)]			
		}else if(type=="seq"){
			file <- f[grep("exp_seq", f)]
		}else{
			stop("Unknown 'type', should be one of: array, seq")
		}
	}
	if(length(file)==0){
		stop("No suitable input files located in the working directory")
	}else if(length(file)>1){
		stop("Multiple suitable file types located in the working directory, please indicate manually using 'file' which file to read and process")
	}
	# Read the input file
	tmp <- read.table(file, sep="\t", header=TRUE)
	if(!verb==0) print(paste("Input file of nrow", nrow(tmp), "read, processing..."))
	patients <- unique(tmp[,"icgc_donor_id"])
	genes <- unique(tmp[,"gene_id"])
	# Create an empty expression matrix and start populating it based on the input data
	m <- matrix(NA, nrow=length(genes), ncol=length(patients))
	colnames(m) <- patients
	rownames(m) <- genes
	# Normalized expression value is platform specific
	expcol <- ifelse(type=="array", 
		# For microarrays use pre-normalized expression
		"normalized_expression_value",
		# For RNA-seq use pre-normalized read counts
		"normalized_read_count"
	)
	# Populate the rows/columns
	if(!verb==0) print(paste("Starting to populate expression matrix of", nrow(m), "genes and", ncol(m), "patients"))
	for(r in 1:nrow(m)){
		if(r %% verb == 0) print(paste("Row number", r, "of", nrow(m)))
		for(c in 1:ncol(m)){
			w <- which(tmp[,"gene_id"] == rownames(m)[r] & tmp[,"icgc_donor_id"] == colnames(m)[c])
			if(length(w)==1){
				m[r,c] <- tmp[w,expcol]
			}else{
				m[r,c] <- mean(tmp[w,expcol])
				## Removed: Do not throw warnings for collapsing multiple probes using mean
				#warning(paste("Multiple instances at ", r, c, collapse=" "))
			}
		}
	}
	m
}

#' Process copy_number_somatic_mutation.*.tsv from ICGC to CNA compatible files; obtain segment mean/median
.mapICGCcn <- function(
	file # A data.frame as read
){
	# If input file/location is not specified try to look for it automatically in current working directory
	if(missing(file)){
		f <- list.files()
		file <- f[grep("copy_number_somatic_mutation", f)]
	}
	if(length(file)==0){
		stop("No suitable input files located in the working directory")
	}else if(length(file)>1){
		stop("Multiple suitable file types located in the working directory, please indicate manually using 'file' which file to read and process")
	}
	# Read the input file
	tmp <- read.table(file, sep="\t", header=TRUE)
	patients <- unique(tmp[,"icgc_donor_id"])
	
}

#' Process copy_number_somatic_mutation.*.tsv from ICGC to CNA compatible files; obtain one of the three originally reported values: "copy neutral LOH" = 0, "loss" = -1, "gain" = +1
#' NOTE:
#' GISTIC-like discretization for the loss/gain appears to differ between datasets;
#' > as.character(unique(PRAD.CA.CNA[,"mutation_type"]))
#' [1] "copy neutral LOH" "loss" "gain"
#' > as.character(unique(PRAD.FR.CNA[,"mutation_type"]))
#' [1] "undetermined"
#' > as.character(unique(PRAD.UK.CNA[,"mutation_type"]))
#' [1] "gain" "amp LOH" "copy neutral LOH" "hemizygous del LOH" "loss" "copy neutral"
#' NOTE: Abbreviation LOH = Loss of heterozygosity
.mapICGCcopynumberGene <- function(
	file # A data.frame as read
){
	# If input file/location is not specified try to look for it automatically in current working directory
	if(missing(file)){
		f <- list.files()
		file <- f[grep("copy_number_somatic_mutation", f)]
	}
	if(length(file)==0){
		stop("No suitable input files located in the working directory")
	}else if(length(file)>1){
		stop("Multiple suitable file types located in the working directory, please indicate manually using 'file' which file to read and process")
	}
	# Read the input file
	tmp <- read.table(file, sep="\t", header=TRUE)
	# Patient IDs
	patients <- unique(tmp[,"icgc_donor_id"])
	tmp <- data.frame(chr = paste("chr", tmp[,"chromosome"], sep=""), start = tmp[,"chromosome_start"], end = tmp[,"chromosome_end"], mutation_type = tmp[,"mutation_type"], icgc_donor_id = tmp[,"icgc_donor_id"])
	# Using GenomicRanges to identify nearest gene
	#if(!require(GenomicFeatures)){
	#	stop("Require 'GenomicFeatures' for GenomicRanges")
	#}
	#hg38_genes <- makeTxDbFromUCSC(genome="hg38", table="refGene")
	#refseq_genes <- genes(hg38_genes)
	## Update, using ACME::findClosestGene instead
	# Using GenomicRanges to identify nearest gene
	if(!require(ACME)){
		stop("Require 'ACME' for findClosestGene() from Bioconductor")
	}
	genenames <- apply(tmp, MARGIN=1, FUN=function(z){
		unique(ACME::findClosestGene(chrom = z["chr"], pos = as.numeric(z["start"]), genome = "hg38", position = "txStart")[,"geneName"])[1]
	})
	tmp[,"gene_id"] <- genenames
	# Create a discretized GISTIC-like CNA matrix with patients on rows and samples as columns
	m <- matrix(NA, nrow=length(genenames), ncol=length(patients))
	rownames(m) <- genenames
	colnames(m) <- patients
	for(r in 1:nrow(m)){
		for(c in 1:ncol(m)){
			w <- which(tmp[,"gene_id"] == rownames(m)[r] & tmp[,"icgc_donor_id"] == colnames(m)[c])
			# Match mutation types and their interpretations for GISTIC-like discretized comparison
			entry <- c(NA, NA, 0, 0, 1, 2, -1, -2)[
				match(tmp[w,"mutation_type"],
					c(
						NA, # = NA
						"undetermined", # = NA
						"copy neutral LOH", # = 0
						"copy neutral", # = 0
						"amp LOH", # = +1
						"gain",  # = +2
						"hemizygous del LOH", # = -1
						"loss" # = -2
					)
				)
			]
			if(length(entry)==1) { m[r,c] <- entry } # Only fill if there is a hit in gene/patient combination in original data
		}
	}
	m
}

###
#
# GISTIC 
# log2 FC <-> GISTIC 2.0 discrete notation for CNA
#
###

# Export rCGH objects to be processed using GISTIC 2.0 in order to harmonize data with cBioPortal CNAs
# 
# From GISTIC 2.0 documentation 
# (URL1: https://cbioportal.readthedocs.io/en/latest/Data-Loading-Tips-and-Best-Practices.html )
# (URL2: ftp://ftp.broadinstitute.org/pub/GISTIC2.0/GISTICDocumentation_standalone.htm )
# (URL3: https://cloud.genepattern.org/gp/pages/login.jsf )
# -> Conflicting information on input files required?
#
# the tab-delimited input segmentation file requires:
#The column headers are: 
#(1)  Sample           (sample name)
#(2)  Chromosome  (chromosome number)
#(3)  Start Position  (segment start position, in bases)
#(4)  End Position   (segment end position, in bases)
#(5)  Num markers      (number of markers in segment)
#(6)  Seg.CN       (log2() -1 of copy number)
.exportGISTIC <- function(
	x, # Should be a list of rCGH-objects for which rCGH:segmentCGH has been run, then segmentation file extracted using getSegTable
	file = "inputGISTIC.tsv" # Output file name for GISTIC 2.0
){
	try({
		if(!class(x)=="list" | !all(lapply(x, FUN=class)=="rCGH-Agilent")){
			stop("Input should be a list of rCGH-objects")
		}
		outputs <- lapply(x, FUN=function(z) {
			tmp <- getSegTable(z)[,c("ID", "chrom", "loc.start", "loc.end", "num.mark", "seg.mean")]			
			tmp[,"ID"] <- z@info["sampleName"]
			tmp
		})
		output <- do.call("rbind", outputs)
		colnames(output) <- c("Sample", "Chrom", "Start", "Stop", "NumMark", "Seg.CN")
		write.table(output, file=file, sep="\t", row.names=F, col.names=T)
	})
}

# Create an integer {-2,-1,0,1,2} matrix for CVA from manually processed GISTIC 2.0 with similar interpretation as processed TCGA data etc
.importGISTIC <- function(
	# GISTIC 2.0 output file as obtained from a GenePattern run; default name here from the result package
	file = "all_thresholded.by_genes.txt"
){
	tmp <- read.table(file, sep="\t", header=TRUE)
	# Omit fields that are no longer required
	tmp <- tmp[,-which(colnames(tmp) %in% c("Locus.ID", "Cytoband"))]
	# Fetch rownames for the new matrix and prune |-suffix
	rnames <- as.character(tmp[,"Gene.Symbol"])
	rnames <- unlist(lapply(rnames, FUN=function(z) { strsplit(z, '|', fixed=TRUE)[[1]][1] }))
	# Cast into an integer matrix with corresponding column names
	tmp <- tmp[,-which(colnames(tmp) == "Gene.Symbol")]
	cnames <- colnames(tmp)
	tmp <- as.matrix(tmp)
	# Integers are sufficient to represent the 5 possible values
	class(tmp) <- "integer"
	# Gene / locus names pruned
	rownames(tmp) <- rnames
	# Column names i.e. GSM-compatible sample names
	colnames(tmp) <- cnames
	# Return the newly formatted matrix
	tmp
}
