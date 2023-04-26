## Steps involved in making the curatedPCaData ExperimentHub-objects

## From the package root:

# ../data-raw/download-clinical.R
## -> Downloads the clinical data matrices from multiple sources (GEO, cBioPortal, cgdsr, ...)
##   and combines them into colData-suitable data.frames for MultiAssayExperiments.

# ../data-raw/download-data.R
## -> Download and harmonize raw data; for example, for gene expression CEL data this will call the 
##    correct functions for using GEOquery, and subsequent data processing such as RMA-normalization.

# ../data-raw/data-curatedPCaData_genes.R
## -> Extract hugo symbols from biomaRt, as well as specific array annotations if available.
##    These are used when constructing e.g. gene-level expression matrices.

# At this stage, the version 1 of curatedPCaData should be build using R CMD build.
# The correct MAE-objects have been stored and in prior versions were offered as LazyLoad'ed MAE-objects.
# After this, derived variables are extracted

# ../data-raw/derive-variables.R
## -> This script assumes the correct MAE-objects are stored inside ../data/-folder of the package.
##    The MAE-objects are loaded from there, and then the derived variables (such as immune deconvolution
##    results and risk scores) are calculated and appended as slots/assays into the existing MAE-objects.

# After above steps, original MultiAssayExperiments were stored in ../data/ of the package.
# In order to make the package compatible with ExperimentHub, these MAE-objects were dissected and an
# additional metadata sheet is automatically generated:

# write.csv(curatedPCaData:::export_metadata(timestamp = "20230215"), file="metadata.csv", quote=TRUE, row.names=FALSE)
## -> Creates the ExperimentHub-friendly metadata.csv

# 