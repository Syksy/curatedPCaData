# tcga data ------
gex_tcga <- curatedPCaData:::generate_cbioportal()
save(gex_tcga, file="data-raw/gex_tcga.RData")

cna_tcga <- curatedPCaData:::generate_cbioportal(
  geneticProfiles="prad_tcga_pub_gistic",
  caseList="prad_tcga_pub_sequenced"
)
save(cna_tcga, file="data-raw/cna_tcga.RData")

# sun et al data -----
gex_sun <- curatedPCaData:::generate_gex_geo()
save(gex_sun, file="data-raw/gex_sun.RData")

#taylor et al data -----
# insert taylor geo data code from daniel here 
# save it to the data-raw folder using commands similar to above


# gex_taylor <- curatedPCaData:::generate_cbioportal(
#   genes = tcga_gene_names,
#   geneticProfiles="prad_mskcc_mrna_median_Zscores",
#   caseList="prad_mskcc_sequenced"
# )
# 
# cna_taylor <- curatedPCaData:::generate_cbioportal(
#   genes = tcga_gene_names,
#   geneticProfiles="prad_mskcc_cna",
#   caseList="prad_mskcc_sequenced"
# )

