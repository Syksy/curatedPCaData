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
# mae_tcga <- curatedPCaData:::create_mae(study_name = "tcga")
# usethis::use_data(mae_tcga, internal = FALSE, overwrite = TRUE)

mae_tcga <- create_mae(study_name = "TCGA")
usethis::use_data(mae_tcga, overwrite = TRUE)

# sun et al data -----
# GEX
gex_sun <- generate_gex_geo(
  geo_code = "GSE25136"
)
save(gex_sun, file="data-raw/gex_sun.RData")

# Create MAE object
# mae_sun <- curatedPCaData:::create_mae(study_name = "sun")
# usethis::use_data(mae_sun, internal = FALSE, overwrite = TRUE)

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
cna_hieronymus <- curatedPCaData:::generate_gex_geo(
  geo_code = "GSE54691"
)
save(cna_hieronymus, file="data-raw/cna_hieronymus.RData")

# Create MAE object

