# tcga data ------
# GEX
gex_tcga <- curatedPCaData:::generate_cbioportal(
  genes = gene_names, # All unique gene symbols
  geneticProfiles = "prad_tcga_pub_rna_seq_v2_mrna", # Omics profile
  caseList = "prad_tcga_pub_sequenced" # Case list
)
usethis::use_data_raw(gex_tcga, overwrite = TRUE)
# CNA
cna_tcga <- curatedPCaData:::generate_cbioportal(
  genes = gene_names,
  geneticProfiles="prad_tcga_pub_gistic",
  caseList="prad_tcga_pub_sequenced"
)
usethis::use_data_raw(cna_tcga, overwrite = TRUE)

# sun et al data -----
# GEX
gex_sun <- curatedPCaData:::generate_gex_geo(
  geo_code = "GSE25136"
)
usethis::use_data_raw(gex_sun, overwrite = TRUE)

#taylor et al data -----
# cBioPortal variant
# GEX
gex_cbio_taylor <- curatedPCaData:::generate_cbioportal(
  genes = gene_names,
  geneticProfiles="prad_mskcc_mrna_median_Zscores",
  caseList="prad_mskcc_sequenced"
)
#usethis::use_data_raw(gex_taylor, overwrite = TRUE)
usethis::use_data_raw(gex_cbio_taylor, overwrite = TRUE)
# CNA
cna_cbio_taylor <- curatedPCaData:::generate_cbioportal(
  genes = gene_names,
  geneticProfiles="prad_mskcc_cna",
  caseList="prad_mskcc_sequenced"
)
#usethis::use_data_raw(cna_taylor, overwrite = TRUE)
usethis::use_data_raw(cna_cbio_taylor, overwrite = TRUE)

# GEO variant
# GEX
gex_geo_taylor <- curatedPCaData:::generate_gex_geo(
  geo_code = "GSE25136"
)
usethis::use_data_raw(gex_geo_taylor, overwrite = TRUE)
# CNA
cna_geo_taylor <- curatedPCaData:::generate_cna_geo(
  geo_code = "GSE21035"
)
usethis::use_data_raw(cna_geo_taylor, overwrite = TRUE)

#hieronymus et al data -----
# CNA
cna_hieronymus <- curatedPCaData::::generate_gex_geo(
  geo_code = "GSE54691"
)
usethis::use_data_raw(cna_hieronymus, overwrite = TRUE)

