###
#
# Sun, et al. portion
#
###

#' Generate gene expression data for Sun et al. 
#' 
#' @param file_directory Character string indicating directory for downloading files. 
#' Files may be large, so please check that there is sufficient free space. 
#' If NULL then files are downloaded into current directory.
#' @param cleanup Logical. Remove tarballs and other files from working directory.
#' @param ... Additional arguments. 
#' @return Gene expression object of a particular type 
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
	keys <- AnnotationDbi::mappedkeys(hgu133a.db::hgu133aGENENAME)
	nam <- names(as.character(hgu133a.db::hgu133aALIAS2PROBE)[match(rownames(GEX_Sun), as.character(hgu133a.db::hgu133aALIAS2PROBE))])
	nam[is.na(nam)] <- "NA"
	# Collapse probes
	GEX_Sun <- do.call("rbind", by(as.matrix(affy::exprs(GEX_Sun)), INDICES=nam, FUN=collapseFUN))
	# Remove downloaded files
	if(cleanup){
		# First GEO download
		file.remove(rownames(supfiles))
		# Tarballs
		file.remove(supfiles2)
	}
	# Return numeric matrix
	as.matrix(GEX_Sun)
}

#' Generate MultiAssaYExperiment-object for Sun, et al.
Generate_MAE_Sun <- function(
	pdata = "../data-raw/GSE25136_curated_pdata.txt",
	image.file,	# If user wishes to save an .RData object of the MAE-object, a file name can be specified
	cleanup,	# Whether intermediate files ought to be cleaned
	...
){
	# Generate gene expression matrix from GEO
	GEX_Sun <- Generate_GEX_Sun(cleanup = cleanup)
	# Start constructing colData
	tmpColDat <- S4Vectors::DataFrame(
		read.table(pdata, # Curated data dictionary file
		sep="\t", header=TRUE)
	)
	rownames(tmpColDat) <- tmpColDat$sample_name
	# Generate a MAE-object for Sun GEX
	MAE_Sun <- MultiAssayExperiment::MultiAssayExperiment(
		experiments = MultiAssayExperiment::ExperimentList(
			"GEX" = GEX_Sun
		),
		colData = tmpColDat,
		sampleMap <- S4Vectors::DataFrame(
			rbind(
				data.frame(assay = "GEX", primary = tmpColDat$patient_id, colname = tmpColDat$sample_name)
			)
		)
	)
	# Optionally save an external R workspace file with the MAE object
	if(!missing(image.file)) save(MAE_Sun, file=image.file)
	# Return MAE object
	MAE_Sun
}


####
#
# Taylor, et al. portion
#
###


#' Gene expression data from Taylor et al.
#'
#' GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE21032
#' Available data types; Affymetrix Exon arrays, 2 separate types
Generate_GEX_Taylor <- function(
	file_directory, 
	cleanup = FALSE,
	...
){
	if(!missing(file_directory)) setwd(file_directory) # exchange setwd with here::here()
	# Fetch supplementary files from Taylor et al. from GEO
	# Be wary as the tarball already is 24.2 Gb
	supfiles <- GEOquery::getGEOSuppFiles('GSE21032')
	# Step in the internal directory created by GEOquery
	setwd("GSE21032")
	# Open the tarball
	utils::untar(tarfile=rownames(supfiles))
	# First download CEL files from GEO
	CELs <- read.celfiles(list.celfiles())	
	#> CELs
	#ExonFeatureSet (storageMode: lockedEnvironment)
	#assayData: 6553600 features, 370 samples 
	#  element names: exprs 
	#protocolData
	#  rowNames: GSM526134_YX_Exon1_PCA0001.CEL GSM526135_YX_Exon1_PCA0002.CEL ... GSM528049_YX_Exon1_PAN0174.CEL (370
	#    total)
	#  varLabels: exprs dates
	#  varMetadata: labelDescription channel
	#phenoData
	#  rowNames: GSM526134_YX_Exon1_PCA0001.CEL GSM526135_YX_Exon1_PCA0002.CEL ... GSM528049_YX_Exon1_PAN0174.CEL (370
	#    total)
	#  varLabels: index
	#  varMetadata: labelDescription channel
	#featureData: none
	#experimentData: use 'experimentData(object)'
	#Annotation: pd.huex.1.0.st.v2
	RMAs <- rma(CELs)
	featureData(RMAs) <- getNetAffx(RMAs, "transcript")
	sampleNames(RMAs) <- unlist(lapply(strsplit(list.celfiles(), "_"), FUN=function(z) z[[1]])) # GSM######-type names from GEO
	nameMapTaylor <- data.frame(
		GSM = unlist(lapply(strsplit(list.celfiles(), "_"), FUN=function(z) z[[1]])), # GEO, GSM-names
		ID = gsub(".CEL", "", unlist(lapply(strsplit(list.celfiles(), "_"), FUN=function(z) z[[4]]))) # mapping to PCA0001, PCA0002, ...
	)
	# 
	GEX_Taylor <- RMAs
	# Proper naming
	nam <- fData(GEX_Taylor)[,8]
	nam2 <- unlist(lapply(nam, FUN=function(z) { strsplit(z, " // ")[[1]][2] }))
	nam2[is.na(nam2)] <- "NA"
	rownames(GEX_Taylor) <- make.unique(nam2)
	## TODO:
	# cleanup
	
	# Return constructed gene expression matrix
	GEX_Taylor
}

#' Copy number alterations from Taylor et al.
#'
#' GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE21032
#' Agilent aCGH
#' Platform in GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL4091
#' Agilent-014693 Human Genome CGH Microarray 244A
#' Ref: Taylor BS, Schultz N, Hieronymus H, Gopalan A et al. Integrative genomic profiling of human prostate cancer. Cancer Cell 2010 Jul 13;18(1):11-22. 
Generate_CNA_Taylor <- function(
	file_directory, 
	cleanup = FALSE,
	...
){
	if(!missing(file_directory)) setwd(file_directory) # exchange setwd with here::here()
	# Fetch supplementary files from Taylor et al. CNA from GEO
	# Be wary as the tarball already is 14.9 Gb
	supfiles <- GEOquery::getGEOSuppFiles('GSE21035')
	# Step in the internal directory created by GEOquery
	setwd("GSE21035")
	# Open the tarball
	utils::untar(tarfile=rownames(supfiles))
	# Read the raw Agilent CGH data
	# Name samples according to the filename
	# Agilent-014693 Human Genome CGH Microarray 244A (Feature number version)
	# rCGH says for readAgilent-function: "Agilent from 44 to 400K are supported."
	TaylorCGH <- lapply(list.files(), FUN=function(z) { 
		try({
			cat("\n\nProcessing: ",z,"\n\n"); 
			rCGH::readAgilent(z, genome="hg38", sampleName=gsub(".txt", "", z)) 
		})
	})
	#> list.files()[which(unlist(lapply(TaylorCGH, FUN=class))=="try-error")]
	#[1] "GSM525755.txt" "GSM525763.txt"
	# Some files appear broken; missing columns?
	
	# Signal adjustments
	TaylorCGH <- lapply(TaylorCGH, FUN=function(z){
		try({
			rCGH::adjustSignal(z) 
		})
	})
	# Segmentation
	TaylorCGH <- lapply(TaylorCGH, FUN=function(z){
		try({
			rCGH::segmentCGH(z) 
		})
	})
	# EM-algorithm normalization
	TaylorCGH <- lapply(TaylorCGH, FUN=function(z){
		try({
			rCGH::EMnormalize(z) 
		})
	})
	# For the old pipeline with lacking sampleName but included fileName, run:
	CNA_Taylor <- lapply(CNA_Taylor, FUN=.renameSampleName)
	### TODO:
	
	# cleanup
	
	## Export a compiled GISTIC 2.0 -compatible file from the segmented samples
	##exportGISTIC(CNA_Taylor, file="inputGISTIC_Taylor.tsv")
	CNA_Taylor
	
}	

#' Generate MultiAssaYExperiment-object for Taylor, et al.
Generate_MAE_Taylor <- function(
	pdata = "../data-raw/prad_mskcc_curated_pdata.txt",
	image.file,	# If user wishes to save an .RData object of the MAE-object, a file name can be specified
	cleanup,	# Whether intermediate files ought to be cleaned
	...
){
	# Gene expression from GEO
	GEX_Taylor <- Generate_GEX_Taylor(cleanup=cleanup)
	# Copy number alterations from GEO
	CNA_Taylor <- Generate_CNA_Taylor(cleanup=cleanup)
	
	# Running Taylor MAE:
	
	# Start constructing colData
	tmpColDat <- S4Vectors::DataFrame(
		#read.table(pdata, # Curated data dictionary file
		read.table("D:/Gits/curatedPCaData/data-raw/prad_mskcc_curated_pdata.txt", # Curated data dictionary file
		sep="\t", header=TRUE)
	)
	# Set row names
	rownames(tmpColDat) <- tmpColDat$patient_id

	# Generate a MAE-object for TCGA 'omics
	MAE_Taylor <- MultiAssayExperiment::MultiAssayExperiment(
		experiments = MultiAssayExperiment::ExperimentList(
			"GEX" = t(GEX_Taylor)#,
			#"CNA" = t(CNA_Taylor)
		),
		colData = tmpColDat,
		sampleMap <- S4Vectors::DataFrame(
			rbind(
				data.frame(assay = "GEX", primary = tmpColDat$patient_id, colname = tmpColDat$sample_name)#,
				#data.frame(assay = "CNA", primary = tmpColDat$patient_id, colname = tmpColDat$sample_name)
			)
		)
	)	
	
	# Optionally save an external R workspace file with the MAE object
	if(!missing(image.file)) save(MAE_Taylor, file=image.file)
	# Return MAE-object from Taylor et al.
	MAE_Taylor
}

####
#
# TCGA portion
#
###

#' Generate a MultiAssayExperiment-object of TCGA GEX + CNA
Generate_MAE_TCGA <- function(
	pdata = "../data-raw/prad_tcga_curated_pdata.txt",	# Location of the curated pdata txt-file
	image.file,	# If user wishes to save an .RData object of the MAE-object, a file name can be specified
	# clean,	# No clean-up required for TCGA as no intermediate files are generated
	...
){
	# Running TCGA GEX:
	GEX_TCGA <- Generate_cBioPortal(genes = genes$hgnc, geneticProfiles="prad_tcga_rna_seq_v2_mrna", caseList="prad_tcga_sequenced")

	# Running TCGA CNA:
	CNA_TCGA <- Generate_cBioPortal(genes = genes$hgnc, geneticProfiles="prad_tcga_gistic", caseList="prad_tcga_sequenced")

	# Running TCGA MAE:
	
	# Start constructing colData
	tmpColDat <- S4Vectors::DataFrame(
		read.table(pdata, # Curated data dictionary file
		sep="\t", header=TRUE)
	)
	# Sanitize '-'-symbols in TCGA names as R tends to replace mathematical operators or whitespace with '.'s
	tmpColDat$sample_name <- gsub("-", ".", tmpColDat$sample_name, fixed=TRUE)
	# Set row names
	rownames(tmpColDat) <- tmpColDat$patient_id

	# Generate a MAE-object for TCGA 'omics
	MAE_TCGA <- MultiAssayExperiment::MultiAssayExperiment(
		experiments = MultiAssayExperiment::ExperimentList(
			"GEX" = t(GEX_TCGA),
			"CNA" = t(CNA_TCGA)
		),
		colData = tmpColDat,
		sampleMap <- S4Vectors::DataFrame(
			rbind(
				data.frame(assay = "GEX", primary = tmpColDat$patient_id, colname = tmpColDat$sample_name),
				data.frame(assay = "CNA", primary = tmpColDat$patient_id, colname = tmpColDat$sample_name)
			)
		)
	)
	# Optionally save an external R workspace file with the MAE object
	if(!missing(image.file)) save(MAE_TCGA, file=image.file)
	# Return MAE object
	MAE_TCGA
}

#' Generic function for pulling data from cBioPortal; useful for e.g. TCGA pulling on multiple omics layers
Generate_cBioPortal <- function(
	genes, # List of gene symbols to iterate over
	geneticProfiles="prad_tcga_rna_seq_v2_mrna", # for cgdsr calls, platform and dataset specific string
	caseList="prad_tcga_sequenced", # for cgdsr calls, platform and dataset specific string
	...
){
	#http://www.cbioportal.org/study?id=prad_tcga#summary
	mycgds <- cgdsr::CGDS("http://www.cbioportal.org/")
	# Use the wrapper function to iterative calls for split gene lists
	dat <- getProfileDataWrapper(
		x=mycgds, # cgdsr object
		#genes=curatedPCaData:::.getGeneNames()$hgnc, # All unique gene symbols
		genes=genes, # All unique gene symbols
		geneticProfiles=geneticProfiles, # Omics profile
		caseList=caseList # Case list
	)
	# TODO: Filter out low quality samples based on list provided by Travis
	# Return omics matrix
	dat
}

###
#
# ICGC portion
#
###

## Multiple specific data pulls from the ICGC depo

#' PRAD-CA (PRAD-CA Prostate Adenocarcinoma - CA)
Generate_ICGC_CA <- function(
	file_directory, 
	# Relevant available files for the Canadian ICGC dataset
	filelist = c(
		"https://dcc.icgc.org/api/v1/download?fn=/release_27/Projects/PRAD-CA/copy_number_somatic_mutation.PRAD-CA.tsv.gz",
		"https://dcc.icgc.org/api/v1/download?fn=/release_27/Projects/PRAD-CA/donor.PRAD-CA.tsv.gz",
		"https://dcc.icgc.org/api/v1/download?fn=/release_27/Projects/PRAD-CA/donor_exposure.PRAD-CA.tsv.gz",
		"https://dcc.icgc.org/api/v1/download?fn=/release_27/Projects/PRAD-CA/donor_family.PRAD-CA.tsv.gz",
		"https://dcc.icgc.org/api/v1/download?fn=/release_27/Projects/PRAD-CA/donor_therapy.PRAD-CA.tsv.gz",
		"https://dcc.icgc.org/api/v1/download?fn=/release_27/Projects/PRAD-CA/exp_array.PRAD-CA.tsv.gz",
		"https://dcc.icgc.org/api/v1/download?fn=/release_27/Projects/PRAD-CA/meth_array.PRAD-CA.tsv.gz",
		"https://dcc.icgc.org/api/v1/download?fn=/release_27/Projects/PRAD-CA/sample.PRAD-CA.tsv.gz",
		"https://dcc.icgc.org/api/v1/download?fn=/release_27/Projects/PRAD-CA/simple_somatic_mutation.open.PRAD-CA.tsv.gz",
		"https://dcc.icgc.org/api/v1/download?fn=/release_27/Projects/PRAD-CA/specimen.PRAD-CA.tsv.gz",
		"https://dcc.icgc.org/api/v1/download?fn=/release_27/Projects/PRAD-CA/structural_somatic_mutation.PRAD-CA.tsv.gz"
	)
){
	if(!missing(file_directory)) setwd(file_directory) # exchange setwd with here::here()
	ICGC.PRAD.CA <- lapply(filelist, FUN=.icgcDownload)
	ICGC.PRAD.CA
}

## PRAD-CN (PRAD-CN Prostate Cancer - CN)
## -> OMIT?

#' PRAD-FR (PRAD-FR Prostate Cancer - Adenocarcinoma - FR)
Generate_ICGC_FR <- function(
	file_directory, 
	# Relevant available files for the Canadian ICGC dataset
	filelist = c(
		"https://dcc.icgc.org/api/v1/download?fn=/release_27/Projects/PRAD-FR/copy_number_somatic_mutation.PRAD-FR.tsv.gz",
		"https://dcc.icgc.org/api/v1/download?fn=/release_27/Projects/PRAD-FR/donor.PRAD-FR.tsv.gz",
		"https://dcc.icgc.org/api/v1/download?fn=/release_27/Projects/PRAD-FR/donor_family.PRAD-FR.tsv.gz",
		"https://dcc.icgc.org/api/v1/download?fn=/release_27/Projects/PRAD-FR/donor_surgery.PRAD-FR.tsv.gz",
		"https://dcc.icgc.org/api/v1/download?fn=/release_27/Projects/PRAD-FR/exp_array.PRAD-FR.tsv.gz",
		"https://dcc.icgc.org/api/v1/download?fn=/release_27/Projects/PRAD-FR/exp_seq.PRAD-FR.tsv.gz",
		"https://dcc.icgc.org/api/v1/download?fn=/release_27/Projects/PRAD-FR/sample.PRAD-FR.tsv.gz",
		"https://dcc.icgc.org/api/v1/download?fn=/release_27/Projects/PRAD-FR/simple_somatic_mutation.open.PRAD-FR.tsv.gz",
		"https://dcc.icgc.org/api/v1/download?fn=/release_27/Projects/PRAD-FR/specimen.PRAD-FR.tsv.gz",
		"https://dcc.icgc.org/api/v1/download?fn=/release_27/Projects/PRAD-FR/structural_somatic_mutation.PRAD-FR.tsv.gz"
	)
){
	if(!missing(file_directory)) setwd(file_directory) # exchange setwd with here::here()
	ICGC.PRAD.FR <- lapply(filelist, FUN=.icgcDownload)
	ICGC.PRAD.FR
}


#' PRAD-UK (PRAD-UK Prostate Adenocarcinoma - UK)
Generate_ICGC_UK <- function(
	file_directory, 
	# Relevant available files for the Canadian ICGC dataset
	filelist = c(
		"https://dcc.icgc.org/api/v1/download?fn=/release_27/Projects/PRAD-UK/copy_number_somatic_mutation.PRAD-UK.tsv.gz",
		"https://dcc.icgc.org/api/v1/download?fn=/release_27/Projects/PRAD-UK/donor.PRAD-UK.tsv.gz",
		"https://dcc.icgc.org/api/v1/download?fn=/release_27/Projects/PRAD-UK/donor_exposure.PRAD-UK.tsv.gz",
		"https://dcc.icgc.org/api/v1/download?fn=/release_27/Projects/PRAD-UK/donor_family.PRAD-UK.tsv.gz",
		"https://dcc.icgc.org/api/v1/download?fn=/release_27/Projects/PRAD-UK/donor_therapy.PRAD-UK.tsv.gz",
		"https://dcc.icgc.org/api/v1/download?fn=/release_27/Projects/PRAD-UK/sample.PRAD-UK.tsv.gz",
		"https://dcc.icgc.org/api/v1/download?fn=/release_27/Projects/PRAD-UK/simple_somatic_mutation.open.PRAD-UK.tsv.gz",
		"https://dcc.icgc.org/api/v1/download?fn=/release_27/Projects/PRAD-UK/specimen.PRAD-UK.tsv.gz",
		"https://dcc.icgc.org/api/v1/download?fn=/release_27/Projects/PRAD-UK/structural_somatic_mutation.PRAD-UK.tsv.gz"
	)
){
	if(!missing(file_directory)) setwd(file_directory) # exchange setwd with here::here()
	ICGC.PRAD.UK <- lapply(filelist, FUN=.icgcDownload)
	ICGC.PRAD.UK
}

## PRAD-US (PRAD-US Prostate Adenocarcinoma - TCGA, US)
## -> OMIT?

####
#
# Supporting variables
#
####

# Extract whole list of desired gene names in various formats
genes <- .getGeneNames()

####
#
# Sun, et al.
# GEX
#
####

# Running Sun, et al.:
MAE_Sun <- Generate_MAE_Sun()

####
#
# TCGA 
# GEX + CNA
#
####

## Get the summary lists for 'omics and case profiles (gives names for omics 'geneticProfiles' and 'caseList')
# cgdsr::getCaseLists(mycgds,"prad_tcga")
# cgdsr::getGeneticProfiles(mycgds,"prad_tcga")

# Note from Jordan & Christelle 
# N=333 in TCGA provisional, Cell 2015

MAE_TCGA <- Generate_MAE_TCGA()


####
#
# Taylor, et al.
# GEX + CNA
#
####

