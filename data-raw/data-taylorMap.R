##
#
# Quick and dirty internal .RData for mapping Taylor IDs between the GEX GEO (GSM##), CNA GEO (GSM##) and cBio compatible IDs (PCA##)
# -- TEMPORARY --
#
##

# Raw Taylor data map as present in Samples under:
#
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE21035 (exon/transcript)
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE21034 (aCGH)
#
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE21032 appears broken (probably because of truncation?)
#
rawTaylorMap <- read.table("data-raw/taylorMap.tsv", sep="\t", header=FALSE)
# Leave out cell lines, xenografts etc
rawTaylorMap <- rawTaylorMap[grep("PCA", rawTaylorMap[,2]),]
# Swap column order and name GSM and PCA columns for MAE as 'colname' and 'primary'
rawTaylorMap <- rawTaylorMap[,1:2]
colnames(rawTaylorMap) <- c("colname", "primary")
# Add an assay column
rawTaylorMap <- cbind(assay = NA, rawTaylorMap)
# Assign corresponding assays
rawTaylorMap[grep("aCGH", rawTaylorMap[,"primary"]),"assay"] <- "cna"
rawTaylorMap[grep("exon|transcript", rawTaylorMap[,"primary"]),"assay"] <- "gex"
# Clean up primary names
rawTaylorMap[,"primary"] <- substr(gsub("Prostate tumor | transcript| exon", "", rawTaylorMap[,"primary"]), start=1, stop=7)

# Save cleaned up map for Taylor et al. 'omics
map_taylor <- rawTaylorMap
save(map_taylor, file="data-raw/map_taylor.RData")
