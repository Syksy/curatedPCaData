####
#
# Helper functions for downloading/updating/fetching packages utilized by datasets used for curatedPCaData
# ideally fetched from CRAN, some from Bioconductor, while rest may only be available on e.g. GitHub
#
####


# Install and update packages for curatedPCaData processing
curatedPCaDataPackages <- function(
	verb = 1, 	# Level of verbosity; <1 omits any output
	update = TRUE, 	# Whether BioConductor packages should be updated (encouraged, dependencies may be out of date)
	ask = FALSE, 	# Whether user should be asked to update Bioconductor packages manually
	biocV = "3.10"	# Assumed BiocManager/Bioconductor version
	
	
){
	noteloadBioc <- function(pckgName="") {
		# Optional verbosity, by default printing out brief cat
		if(verb>=1){ cat(paste("Loading/installing '",pckgName,"'...\n"))}
		# Always make sure BiocManager is functional before going to Bioconductor-packages
		if (!requireNamespace("BiocManager", quietly = TRUE))
		{
			install.packages("BiocManager")
			BiocManager::install(version = biocV, update=update, ask=ask)
			# May require manual input depending on parameters provided to 'install':
			#> BiocManager::install(version = '3.10')
			#Upgrade 83 packages to Bioconductor version '3.10'? [y/n]:
			#library("BiocManager")
		}
		# Load/install package of interest
		if(!require(pckgName, character.only=TRUE)){
			BiocManager::install(as.character(pckgName), version = biocV, update=update, ask=ask) 
			library(pckgName, character.only=TRUE)
		}
	}

	##
	#
	# Utilized by multiple pipelines
	#
	##
	
	# For fetching data directly from GEO
	noteloadBioc("GEOquery")

	# Ref: Commo F, Guinney J, Ferte C, Bot B, Lefebvre C, Soria JC, and Andre F.
	# rcgh : a comprehensive array-based genomic profile platform for precision
	# medicine. Bioinformatics, 2015.
	# https://bioconductor.org/packages/release/bioc/vignettes/rCGH/inst/doc/rCGH.pdf
	noteloadBioc("rCGH")

	# GenomicRanges
	noteloadBioc("GenomicFeatures")

	# Unique gene name lists
	noteloadBioc("TxDb.Hsapiens.UCSC.hg38.knownGene")

	# Required by Taylor, et al. (CNA)
	noteloadBioc("BiocParallel")

	# Required by Taylor, et al. (GEX)
	#BiocManager::install("frma", version = "3.8") # No longer in use for the pipeline
	noteloadBioc("affy")
	noteloadBioc("oligo")
}

curatedPCaDataPackages()

