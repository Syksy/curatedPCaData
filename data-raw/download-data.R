# tcga data ------
# GEX
gex_tcga <- generate_cbioportal(
  genes = sort(unique(curatedPCaData_genes$hgnc_symbol)), # All unique gene symbols
  geneticProfiles = "prad_tcga_pub_rna_seq_v2_mrna", # Omics profile
  caseList = "prad_tcga_pub_sequenced" # Case list
)
save(gex_tcga, file="data-raw/gex_tcga.RData")

# CNA
cna_tcga <- generate_cbioportal(
  genes = sort(unique(curatedPCaData_genes$hgnc_symbol)),
  geneticProfiles="prad_tcga_pub_linear_CNA", # changed from GISTIC to linear values to be comparable to log2 FCs from other datasets
  caseList="prad_tcga_pub_sequenced"
)
save(cna_tcga, file="data-raw/cna_tcga.RData")

# Create MAE object
mae_tcga <- create_mae(study_name = "TCGA")
usethis::use_data(mae_tcga, overwrite = TRUE)

# sun et al data -----
# GEX
gex_sun <- generate_gex_geo(
  geo_code = "GSE25136"
)
save(gex_sun, file="data-raw/gex_sun.RData")

# Create MAE object
mae_sun <- create_mae(study_name = "Sun")
usethis::use_data(mae_sun, overwrite = TRUE)

#taylor et al data -----
gex_taylor <- curatedPCaData:::generate_gex_geo(
  geo_code = "GSE25136"
)
save(gex_taylor, file="data-raw/gex_taylor.RData")

# CNA
cna_taylor <- curatedPCaData:::generate_cna_geo(
  geo_code = "GSE21035"
)
save(cna_taylor, file="data-raw/cna_taylor.RData")

# Create MAE object
mae_taylor <- create_mae(study_name = "Taylor")
usethis::use_data(mae_taylor, internal = FALSE, overwrite = TRUE)

#hieronymus et al data -----
# CNA
cna_hieronymus <- curatedPCaData:::generate_cna_geo(
  geo_code = "GSE54691"
)
save(cna_hieronymus, file="data-raw/cna_hieronymus.RData")

# Create MAE object
mae_hieronymus <- curatedPCaData:::create_mae(study_name = "hieronymus")
usethis::use_data(mae_hieronymus, internal = FALSE, overwrite = TRUE)


##
#
# ICGC datasets
#
##

# PRAD-CA
gex_icgcca <- curatedPCaData:::generate_icgc("PRAD_CA", "gex")
save(gex_icgcca, file="data-raw/gex_icgcca.RData")

# Create MAE object
mae_icgcca <- curatedPCaData:::create_mae(study_name = "icgcca")
usethis::use_data(mae_icgcca, internal = FALSE, overwrite = TRUE)

# PRAD-FR
gex_icgcfr <- curatedPCaData::generate_icgc("PRAD_FR", "gex")
save(gex_icgcfr, file="data-raw/gex_icgcfr.RData")

# TODO: At this point the data contains raw read counts, and is not yet usable as a 2-dim gex matrix

# TODO: MAE

# PRAD-UK

###################################################
# Friedrich 2020 FOR NOW BASED ON PROCESSED DATA!!!
###################################################

library(GEOquery)
library(Biobase)

# load series and platform data from GEO

fr_gset <- getGEO("GSE134051", GSEMatrix =TRUE, getGPL=TRUE)

labels = Biobase::fData(fr_gset[[1]])
gtab = curatedPCaData:::curatedPCaData_genes

if (length(fr_gset) > 1) idx <- grep("GPL26898", attr(gset, "names")) else idx <- 1
fr_ex <- exprs(fr_gset[[idx]])

# replacing row names with gene ids
##############################################
labels$ensb = substr(labels$SPOT_ID, 1, 15)
rownames(fr_ex) = labels$ensb
fr_ex = fr_ex[rownames(fr_ex) != 'NoEntry', ]
fr_ex = fr_ex[substr(rownames(fr_ex), 1, 4) != 'XLOC', ]
fr_ex = fr_ex[is.element(rownames(fr_ex), gtab[,1]), ]
gtab2 = gtab[match(rownames(fr_ex), gtab[,1]), ]

gtab2[which(gtab2[,3] == ''), 3] = gtab2[which(gtab2[,3] == ''), 1]

rownames(fr_ex) = gtab2[,3]
gex_friedrich = aggregate(fr_ex, by = list(rownames(ex)), mean)
rownames(gex_friedrich) = gex_friedrich[,1]
gex_friedrich = gex_friedrich[, -c(1)]

save(gex_friedrich, file = "data-raw/gex_friedrich.RData")

mae_friedrich = curatedPCaData:::create_mae(study_name = 'Friedrich')








