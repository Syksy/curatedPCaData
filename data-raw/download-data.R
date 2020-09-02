# tcga data ------
gex_tcga <- curatedPCaData:::generate_cbioportal(
  genes = tcga_gene_names, # All unique gene symbols
  geneticProfiles = "prad_tcga_pub_rna_seq_v2_mrna", # Omics profile
  caseList = "prad_tcga_pub_sequenced" # Case list
)
usethis::use_data_raw(gex_tcga, overwrite = TRUE)

cna_tcga <- curatedPCaData:::generate_cbioportal(
  genes = tcga_gene_names,
  geneticProfiles="prad_tcga_pub_gistic",
  caseList="prad_tcga_pub_sequenced"
)
usethis::use_data_raw(cna_tcga, overwrite = TRUE)

# sun et al data -----
gex_sun <- curatedPCaData:::generate_gex_geo()
usethis::use_data_raw(gex_sun, overwrite = TRUE)

#taylor et al data -----
gex_taylor <- curatedPCaData:::generate_cbioportal(
  genes = tcga_gene_names,
  geneticProfiles="prad_mskcc_mrna_median_Zscores",
  caseList="prad_mskcc_sequenced"
)
usethis::use_data_raw(gex_taylor, overwrite = TRUE)

cna_taylor <- curatedPCaData:::generate_cbioportal(
  genes = tcga_gene_names,
  geneticProfiles="prad_mskcc_cna",
  caseList="prad_mskcc_sequenced"
)
usethis::use_data_raw(cna_taylor, overwrite = TRUE)
