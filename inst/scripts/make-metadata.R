## Steps involved in creating the metadata.csv provided for ExperimentHub

## The code for generating metadata.csv for ExperimentHub are based on the MAE-objects previously present in curatedPCaData; the function call for generating this is:
## Located in 'curatedPCaData/R/exportfuncs.R'
## Example export with time stamp for 2023, February 15th:
write.csv(curatedPCaData:::export_metadata(timestamp = "20230215"), file="metadata.csv", quote=TRUE, row.names=FALSE)

