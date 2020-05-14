####
#
# Helper functions for downloading/updating/fetching packages utilized by datasets used for curatedPCaData
# ideally fetched from CRAN, some from Bioconductor, while rest may only be available on e.g. GitHub
#
####


# Install and update packages for curatedPCaData processing
curatedPCaDataPackages <- function(
	verb = 1, 	# Level of verbosity; <1 omits any output
	biocV = "3.10"	# Assumed BiocManager/Bioconductor version
){
	noteloadBioc <- function(
		pckgName="",	# Name of the package as character string
		update = TRUE, 	# Whether BioConductor packages should be updated (encouraged, dependencies may be out of date)
		ask = FALSE 	# Whether user should be asked to update Bioconductor packages manually
	) {
		# Optional verbosity, by default printing out brief cat
		if(verb>=1){ cat(paste("Loading/installing '",pckgName,"'...\n"))}
		# Always make sure BiocManager is functional before going to Bioconductor-packages
		if (!requireNamespace("BiocManager", quietly = TRUE))
		{
			install.packages("BiocManager")
			BiocManager::install(version = biocV, update = update, ask = ask)
			# May require manual input depending on parameters provided to 'install':
			#> BiocManager::install(version = '3.10')
			#Upgrade 83 packages to Bioconductor version '3.10'? [y/n]:
			#library("BiocManager")
		}
		# Load/install package of interest
		if(!require(pckgName, character.only=TRUE)){
			BiocManager::install(as.character(pckgName), version = biocV, update = update, ask = ask) 
			library(pckgName, character.only=TRUE)
		}
	}

	###
	#
	# Double checking that some CRAN R-packages are available
	# TODO: Move to Requires/Imports/Suggests
	# 
	###

	if(!require("RMariaDB")){
		install.packages("RMariaDB")
		library("RMariaDB")
	}

	if(!require("pillar")){
		install.packages("pillar")
		library("pillar")
	}

	##
	#
	# Utilized by multiple pipelines
	#
	##

	# Some base packages that are required downstream and may be in use if not installed early
	try({noteloadBioc("Biobase")})
	try({noteloadBioc("S4Vectors")})
	try({noteloadBioc("GenomicFeatures")})
	
	# For fetching data directly from GEO
	try({noteloadBioc("GEOquery")})

	# Unique gene name lists / required by rCGH (some may be redundant)
	try({noteloadBioc("TxDb.Hsapiens.UCSC.hg18.knownGene", update = FALSE)})
	try({noteloadBioc("TxDb.Hsapiens.UCSC.hg19.knownGene", update = FALSE)})
	try({noteloadBioc("TxDb.Hsapiens.UCSC.hg38.knownGene", update = FALSE)})

	## Needed for CNA processing
	# Ref: Commo F, Guinney J, Ferte C, Bot B, Lefebvre C, Soria JC, and Andre F.
	# rcgh : a comprehensive array-based genomic profile platform for precision
	# medicine. Bioinformatics, 2015.
	# https://bioconductor.org/packages/release/bioc/vignettes/rCGH/inst/doc/rCGH.pdf
	try({noteloadBioc("rCGH")})

	# GenomicRanges
	try({noteloadBioc("GenomicFeatures")})

	# Required by Taylor, et al. (CNA)
	try({noteloadBioc("BiocParallel")})

	# Required by Taylor, et al. (GEX)
	try({noteloadBioc("oligo")})

	# Taylor, et al.; Sun, et al.
	try({noteloadBioc("affy")}) 
	
	# Required by Sun, et al. (GEX)
	try({noteloadBioc("hgu133a.db")})
	
	# Required by Taylor, et al. (GEX)
	try({noteloadBioc("pd.huex.1.0.st.v2")})
	
	###
	#
	# Generalized, non-dataset or 'omics specific toolkits
	#
	###
	
	# MultiAssayExperiment-objects
	try({noteloadBioc("MultiAssayExperiment")
	
}
# Runnable as:
# > curatedPCaDataPackages()

