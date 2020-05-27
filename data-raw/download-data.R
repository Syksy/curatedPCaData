tcga_gex <- curatedPCaData:::generate_cbioportal(
  genes = tcga_gene_names, # All unique gene symbols
  geneticProfiles = "prad_tcga_pub_rna_seq_v2_mrna", # Omics profile
  caseList = "prad_tcga_pub_sequenced" # Case list
)
usethis::use_data(tcga_gex, overwrite = TRUE)

tcga_cna <- curatedPCaData:::generate_cbioportal(
  genes = tcga_gene_names,
  geneticProfiles="prad_tcga_pub_gistic",
  caseList="prad_tcga_pub_sequenced"
)
usethis::use_data(tcga_cna, overwrite = TRUE)

sun_gex <- curatedPCaData:::generate_gex_geo()
usethis::use_data(sun_gex, overwrite = TRUE)

taylor_gex <- curatedPCaData:::generate_cbioportal(
  genes = tcga_gene_names,
  geneticProfiles="prad_mskcc_mrna_median_Zscores",
  caseList="prad_mskcc_sequenced"
)
usethis::use_data(tcga_gex, overwrite = TRUE)

taylor_cna <- curatedPCaData:::generate_cbioportal(
  genes = tcga_gene_names,
  geneticProfiles="prad_mskcc_cna",
  caseList="prad_mskcc_sequenced"
)
usethis::use_data(tcga_cna, overwrite = TRUE)
