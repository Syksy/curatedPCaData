# tcga data ------
# GEX
gex_tcga <- curatedPCaData:::generate_cbioportal(
  genes = sort(unique(curatedPCaData:::curatedPCaData_genes$hgnc_symbol)), # All unique gene symbols
  geneticProfiles = "prad_tcga_pub_rna_seq_v2_mrna", # Omics profile
  caseList = "prad_tcga_pub_sequenced" # Case list
)
#save(gex_tcga, file="data-raw/gex_tcga.RData")
usethis::use_data_raw(gex_tcga, overwrite = TRUE)
# CNA
cna_tcga <- curatedPCaData:::generate_cbioportal(
  genes = sort(unique(curatedPCaData:::curatedPCaData_genes$hgnc_symbol)),
  geneticProfiles="prad_tcga_pub_linear_CNA", # changed from GISTIC to linear values to be comparable to log2 FCs from other datasets
  caseList="prad_tcga_pub_sequenced"
)
#save(cna_tcga, file="data-raw/cna_tcga.RData")
usethis::use_data_raw(cna_tcga, overwrite = TRUE)
# Create MAE object
mae_tcga <- curatedPCaData:::create_mae(study_name = "tcga")
usethis::use_data(mae_tcga, internal = FALSE, overwrite = TRUE)

# sun et al data -----
# GEX
gex_sun <- curatedPCaData:::generate_gex_geo(
  geo_code = "GSE25136"
)
usethis::use_data_raw(gex_sun, overwrite = TRUE)
#save(gex_sun, file="data-raw/gex_sun.RData")

#taylor et al data -----
# cBioPortal variant
# GEX
gex_cbio_taylor <- curatedPCaData:::generate_cbioportal(
  genes = sort(unique(curatedPCaData:::curatedPCaData_genes$hgnc_symbol)),
  geneticProfiles="prad_mskcc_mrna_median_Zscores",
  caseList="prad_mskcc_sequenced"
)
#usethis::use_data_raw(gex_taylor, overwrite = TRUE)
usethis::use_data_raw(gex_cbio_taylor, overwrite = TRUE)
#save(gex_cbio_taylor, file="data-raw/gex_cbio_taylor.RData")
# CNA
cna_cbio_taylor <- curatedPCaData:::generate_cbioportal(
  genes = sort(unique(curatedPCaData:::curatedPCaData_genes$hgnc_symbol)),
  geneticProfiles="prad_mskcc_cna",
  caseList="prad_mskcc_sequenced"
)
#usethis::use_data_raw(cna_taylor, overwrite = TRUE)
#save(cna_cbio_taylor, file="data-raw/cna_cbio_taylor.RData")
usethis::use_data_raw(cna_cbio_taylor, overwrite = TRUE)

# GEO variant
# GEX
gex_geo_taylor <- curatedPCaData:::generate_gex_geo(
  geo_code = "GSE25136"
)
#save(gex_geo_taylor, file="data-raw/gex_geo_taylor.RData")
usethis::use_data_raw(gex_geo_taylor, overwrite = TRUE)
# CNA
cna_geo_taylor <- curatedPCaData:::generate_cna_geo(
  geo_code = "GSE21035"
)
#save(cna_geo_taylor, file="data-raw/cna_geo_taylor.RData")
usethis::use_data_raw(cna_geo_taylor, overwrite = TRUE)

#hieronymus et al data -----
# CNA
cna_hieronymus <- curatedPCaData::::generate_gex_geo(
  geo_code = "GSE54691"
)
#save(cna_hieronymus, file="data-raw/cna_hieronymus.RData")
usethis::use_data_raw(cna_hieronymus, overwrite = TRUE)

