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

#Download OSF data
osf_download <- osf_retrieve_file("https://osf.io/m5nh6/") %>% osf_download("./data-raw")
R.utils::gunzip("data-raw/TCGA_PRAD_tpm.tsv.gz")
osf_data <- rio::import("data-raw/TCGA_PRAD_tpm.tsv")

#Download the mapping file 
osf <- curatedPCaData:::format_osf_data(osf_data)
osf_retrieve_file("https://osf.io/7qpsg/")%>% osf_download("./data-raw")

# Re-format the OSF data
osf_t <- t(osf)
osf_t <- as.data.frame(osf_t)
colnames(osf_t) <- osf_t[1,]
osf <- osf_t[-1,]

colnames(osf) <- gsub(x = colnames(osf), pattern = "-", replacement = ".")  
colnames(osf) <- paste(colnames(osf), '01', sep='.')
usethis::use_data(osf, internal = TRUE, overwrite = TRUE)
save(osf, file="data-raw/osfgex_tcga.RData")

# PRAD Barbieri ------
# GEX
gex_barbieri <- generate_cbioportal(
  genes = sort(unique(curatedPCaData_genes$hgnc_symbol)), # All unique gene symbols
  geneticProfiles = "prad_broad_mrna", # Omics profile
  caseList = "prad_broad_sequenced" # Case list
)
save(gex_barbieri, file="data-raw/gex_barbieri.RData")

# CNA
cna_barbieri <- generate_cbioportal(
  genes = sort(unique(curatedPCaData_genes$hgnc_symbol)),
  geneticProfiles="prad_broad_cna", 
  caseList="prad_broad_sequenced"
)
save(cna_barbieri, file="data-raw/cna_barbieri.RData")

# Create MAE object
mae_barbieri <- create_mae(study_name = "Barbieri")
usethis::use_data(mae_barbieri, overwrite = TRUE)

# PRAD Ren ------

# GEX
gex_ren <- generate_cbioportal(
  genes = sort(unique(curatedPCaData_genes$hgnc_symbol)), # All unique gene symbols
  geneticProfiles = "prad_eururol_2017_rna_seq_mrna", # Omics profile
  caseList = "prad_eururol_2017_sequenced" # Case list
)
save(gex_ren, file="data-raw/gex_ren.RData")

# CNA
cna_ren <- generate_cbioportal(
  genes = sort(unique(curatedPCaData_genes$hgnc_symbol)),
  geneticProfiles="prad_eururol_2017_cna", 
  caseList="prad_eururol_2017_sequenced"
)
save(cna_ren, file="data-raw/cna_ren.RData")

# Create MAE object
mae_ren <- create_mae(study_name = "Ren")
usethis::use_data(mae_ren, overwrite = TRUE)



