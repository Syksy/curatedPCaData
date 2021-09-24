### Data downloading scripts for:
## GEX: Gene expression
## CNA: Copy number alteration
## MUT: Mutations
### Available platforms for extracting raw or processed data
## GEO: Gene Omnibus
## cBioPortal (query via cgdsr-package)

### Rough pipeline:
## Download clinical data first (via download-clinical.R), and save in data-raw
## Download GEX/CNA/MUT and save in data-raw using functions via generate.R
## Create final MAE object for package export using functions via create.R


### Alphabetic ordering of datasets:
## - Baca et al.
## - Barbieri et al.
## - Barwick et al. (TODO, not exported yet)
## - Chandran et al.
## - Friedrich et al.
## - Hieronymus et al.
## - ICGC sub datasets
## - IGC
## - Kim et al.
## - Kunderfranco et al.
## - Ren et al.
## - Sun et al.
## - Taylor et al. (MSKCC)
## - TCGA
## - Wallace et al.
## - Wang et al.
## (to be updated)



## - Abida et al. -
# GEX PolyA (FPKM)
gex_polyA_abida <- curatedPCaData:::generate_cbioportal(
  genes = sort(unique(curatedPCaData:::curatedPCaData_genes$hgnc_symbol)), # All unique gene symbols
  geneticProfiles = "prad_su2c_2019_mrna_seq_fpkm_polya_all_sample_Zscores", # Omics profile
  caseList = "prad_su2c_2019_all" # Case list
)
save(gex_polyA_abida, file="data-raw/gex_polyA_abida.RData")

# CNA
cna_abida <- curatedPCaData:::generate_cbioportal(
  genes = sort(unique(curatedPCaData:::curatedPCaData_genes$hgnc_symbol)),
  geneticProfiles="prad_su2c_2019_gistic", 
  caseList="prad_su2c_2019_sequenced"
)
save(cna_abida, file="data-raw/cna_abida.RData")

# Mutations
mut_abida <- curatedPCaData:::generate_cbioportal(
  genes = sort(unique(curatedPCaData:::curatedPCaData_genes$hgnc_symbol)),
  geneticProfiles="prad_su2c_2019_mutations", 
  caseList="prad_su2c_2019_sequenced"
)
mut_abida[which(mut_abida=="NaN")] <- NA
save(mut_abida, file="data-raw/mut_abida.RData")

# Create MAE object
mae_abida <- curatedPCaData:::create_mae(study_name = "abida")
usethis::use_data(mae_abida, overwrite = TRUE)

## - end Abida et al. -


## - Baca et al. -
# CNA
cna_baca <- curatedPCaData:::generate_cbioportal(
  genes = sort(unique(curatedPCaData:::curatedPCaData_genes$hgnc_symbol)),
  geneticProfiles="prad_broad_2013_cna", 
  caseList="prad_broad_2013_sequenced"
)
save(cna_baca, file="data-raw/cna_baca.RData")

# Create MAE object
mae_baca <- create_mae(study_name = "baca")
usethis::use_data(mae_baca, overwrite = TRUE)

## - end Baca et al. -


## - Barbieri et al. -
# GEX
gex_barbieri <- curatedPCaData:::generate_cbioportal(
  genes = sort(unique(curatedPCaData:::curatedPCaData_genes$hgnc_symbol)), # All unique gene symbols
  geneticProfiles = "prad_broad_mrna", # Omics profile
  caseList = "prad_broad_sequenced" # Case list
)
save(gex_barbieri, file="data-raw/gex_barbieri.RData")

# CNA
cna_barbieri <- curatedPCaData:::generate_cbioportal(
  genes = sort(unique(curatedPCaData:::curatedPCaData_genes$hgnc_symbol)),
  geneticProfiles="prad_broad_cna", 
  caseList="prad_broad_sequenced"
)
save(cna_barbieri, file="data-raw/cna_barbieri.RData")

# MUT
mut_barbieri <- curatedPCaData:::generate_cbioportal(
  genes = sort(unique(curatedPCaData:::curatedPCaData_genes$hgnc_symbol)),
  geneticProfiles="prad_broad_mutations", 
  caseList="prad_broad_sequenced"
)
mut_barbieri[which(mut_barbieri=="NaN")] <- NA
save(mut_barbieri, file="data-raw/mut_barbieri.RData")

# Create MAE object
mae_barbieri <- curatedPCaData:::create_mae(study_name = "Barbieri")
usethis::use_data(mae_barbieri, overwrite = TRUE)

## - end Barbieri et al. -


## - Barwick et al.
# GEX: GPL5858	DASL Human Cancer Panel by Gene
gex.logq_barwick <- curatedPCaData:::generate_gex_geo(
	geo_code = "GSE18655"
)
save(gex.logq_barwick, file="data-raw/gex.logq_barwick.RData")

# Create MAE object
mae_barwick <- curatedPCaData:::create_mae(study_name = "Barwick")
usethis::use_data(mae_barwick, overwrite = TRUE)
 
## - end Barwick et al. -


## - Chandran et al. (Yu et al.) -
# GEX: Multiple array types:
# - GPL92	[HG_U95B] Affymetrix Human Genome U95B Array
# - GPL93	[HG_U95C] Affymetrix Human Genome U95C Array
# - GPL8300	[HG_U95Av2] Affymetrix Human Genome U95 Version 2 Array
gex_chandran <- curatedPCaData::generate_gex_geo(
	geo_code = "GSE6919"
)
save(gex_chandran, file="data-raw/gex_chandran.RData")

# Create and save MAE object
mae_chandran <- curatedPCaData:::create_mae(study_name = "chandran")
usethis::use_data(mae_chandran, internal = FALSE, overwrite = TRUE)

# - end Chandran et al. -


## - Friedrich et al. (2020) -
# GEX: GPL26898	Agilent-058029 Custom human expression microarray (Probe Name version)
gex_friedrich <- curatedPCaData:::generate_gex_geo(
	geo_code = "GSE134051"
)
save(gex_friedrich, file="data-raw/gex_friedrich.RData")

mae_friedrich = curatedPCaData:::create_mae(study_name = 'Friedrich')
usethis::use_data(mae_friedrich, internal = FALSE, overwrite = TRUE)

## - end Friedrich et al. -


## - Hieronymus et al. -
# CNA
cna_hieronymus <- curatedPCaData:::generate_cna_geo(
  geo_code = "GSE54691"
)
save(cna_hieronymus, file="data-raw/cna_hieronymus.RData")

# Create MAE object
mae_hieronymus <- curatedPCaData:::create_mae(study_name = "hieronymus")
usethis::use_data(mae_hieronymus, internal = FALSE, overwrite = TRUE)

## - end Hieronymus et al. -


## - ICGC datasets -

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

## - end ICGC datatasets -



## - IGC -
# GEX: [HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array
gex.rma_igc <- curatedPCaData:::generate_gex_geo(
	geo_code = "GSE2109",
	cleanup = FALSE
	# Do not download the batch files, filter the raw dump and filter down to prostate specific samples afterwards
	filter_regex = "_RAW",
)
save(gex.rma_igc, file="data-raw/gex.rma_igc.RData")

# Create MAE object
mae_igc <- curatedPCaData:::create_mae(study_name = "igc")
usethis::use_data(mae_igc, overwrite = TRUE)

## - end IGC -


## - Kim et al. -
# GEX: [HuEx-1_0-st] Affymetrix Human Exon 1.0 ST Array [probe set (exon) version]
gex.rma_kim <- curatedPCaData:::generate_gex_geo(
	geo_code = "GSE119616"
	pckg = "oligo",
	cleanup = FALSE
  
)
save(gex.rma_kim, file="data-raw/gex.rma_kim.RData")

# Create MAE object
mae_kim <- curatedPCaData:::create_mae(study_name = "kim")
usethis::use_data(mae_kim, overwrite = TRUE)

## - end Kim et al. -


## - Kunderfranco et al. -
# GEX: Agilent-012097 Human 1A Microarray (V2) G4110B (Feature Number version)
gex_kunderfranco <- curatedPCaData:::generate_gex_geo(
	geo_code = "GSE14206",
	pckg = "oligo",
	cleanup = FALSE	
)
save(gex_kunderfranco, file = "data-raw/gex_kunderfranco.RData")

# Create MAE object
mae_kunderfranco <- curatedPCaData:::create_mae(study_name = "kunderfranco")
usethis::use_data(mae_kunderfranco, internal = FALSE, overwrite = TRUE)

## - end Kunderfranco et al. -



## - Ren et al. -
# GEX
gex_ren <- curatedPCaData:::generate_cbioportal(
  genes = sort(unique(curatedPCaData:::curatedPCaData_genes$hgnc_symbol)), # All unique gene symbols
  geneticProfiles = "prad_eururol_2017_rna_seq_mrna", # Omics profile
  caseList = "prad_eururol_2017_sequenced" # Case list
)
save(gex_ren, file="data-raw/gex_ren.RData")

# CNA
cna_ren <- curatedPCaData:::generate_cbioportal(
  genes = sort(unique(curatedPCaData:::curatedPCaData_genes$hgnc_symbol)),
  geneticProfiles="prad_eururol_2017_cna", 
  caseList="prad_eururol_2017_sequenced"
)
save(cna_ren, file="data-raw/cna_ren.RData")

# Mutations
mut_ren <- curatedPCaData:::generate_cbioportal(
  genes = sort(unique(curatedPCaData:::curatedPCaData_genes$hgnc_symbol)), # All unique gene symbols
  geneticProfiles = "prad_eururol_2017_mutations", # Omics profile
  caseList = "prad_eururol_2017_sequenced" # Case list
)
mut_ren[which(mut_ren=="NaN")] <- NA
save(mut_ren, file="data-raw/mut_ren.RData")

# Create MAE object
mae_ren <- curatedPCaData:::create_mae(study_name = "ren")
usethis::use_data(mae_ren, overwrite = TRUE)

## - end Ren et al.


## - Sun et al.  -
# GEX: [HG-U133A] Affymetrix Human Genome U133A Array
gex.rma_sun <- curatedPCaData:::generate_gex_geo(
	geo_code = "GSE25136",
	cleanup = FALSE,
	pckg = "oligo"	
)
save(gex.rma_sun, file="data-raw/gex.rma_sun.RData")

# Create MAE object
mae_sun <- curatedPCaData:::create_mae(study_name = "Sun")
usethis::use_data(mae_sun, overwrite = TRUE)

## - end Sun et al. -


## Taylor et al., also known as the MSKCC data -----
# GEX:
# - GPL5188	[HuEx-1_0-st] Affymetrix Human Exon 1.0 ST Array [probe set (exon) version]
# - GPL8227	Agilent-019118 Human miRNA Microarray 2.0 G4470B (miRNA ID version)
# - GPL10264	Affymetrix Human Exon 1.0 ST Array [CDF: HuEx_1_0_st_v2_main_A20071112_EP.cdf]
gex.rma_taylor <- curatedPCaData:::generate_gex_geo(
	geo_code = "GSE21032",
	cleanup = FALSE,
	pckg = "oligo"
)
save(gex_taylor, file="data-raw/gex.rma_taylor.RData")

# CNA
cna_taylor <- curatedPCaData:::generate_cna_geo(
  geo_code = "GSE21035"
)
save(cna_taylor, file="data-raw/cna_taylor.RData")

# Mutations - notice this is downloaded from cBioPortal rather than processed from GEO
mut_taylor <- curatedPCaData:::generate_cbioportal(
  genes = sort(unique(curatedPCaData:::curatedPCaData_genes$hgnc_symbol)),
  geneticProfiles = "prad_mskcc_mutations",
  caseList="prad_mskcc_sequenced"
) 
# Save NA values as truly NA instead of "NaN" even if other instances exist on column
mut_taylor[which(mut_taylor=="NaN")] <- NA
# Grep down to using only patient samples, omitting cell lines etc
mut_taylor <- mut_taylor[,grep("PCA", colnames(mut_taylor))]
save(mut_taylor, file="data-raw/mut_taylor.RData")

# Create MAE object
mae_taylor <- curatedPCaData:::create_mae(study_name = "Taylor")
usethis::use_data(mae_taylor, internal = FALSE, overwrite = TRUE)

## - end Taylor et al. (MSKCC) -


## - TCGA - 
#Download OSF GEX data
osf_download <- osf_retrieve_file("https://osf.io/m5nh6/") %>% osf_download("./data-raw")
R.utils::gunzip("data-raw/TCGA_PRAD_tpm.tsv.gz")
osf_data <- rio::import("data-raw/TCGA_PRAD_tpm.tsv")

#Download the mapping file 
osf <- format_osf_data(osf_data)
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

load("data-raw/osfgex_tcga.RData")
osf<-as.matrix(osf)
storage.mode(osf) <- "numeric"
osf_gex_rounded<-round(osf,digits = 1)
save(osf_gex_rounded, file="data-raw/osfgex_rounded_tcga.RData")
unlink("data-raw/osfgex_tcga.RData")
file.rename("data-raw/osfgex_rounded_tcga.RData","data-raw/gex_tcga.RData")

# GEX
# gex_tcga <- curatedPCaData:::generate_cbioportal(
#   genes = sort(unique(curatedPCaData:::curatedPCaData_genes$hgnc_symbol)), # All unique gene symbols
#   geneticProfiles = "prad_tcga_pub_rna_seq_v2_mrna", # Omics profile
#   caseList = "prad_tcga_pub_sequenced" # Case list
# )
# save(gex_tcga, file="data-raw/gex_tcga.RData")

# CNA
cna_tcga <- curatedPCaData:::generate_cbioportal(
  genes = sort(unique(curatedPCaData:::curatedPCaData_genes$hgnc_symbol)),
  #geneticProfiles="prad_tcga_pub_linear_CNA", # changed from GISTIC to linear values to be comparable to log2 FCs from other datasets
  geneticProfiles="prad_tcga_pub_gistic", # changed back to GISTIC for interpretability in oncoprints
  caseList="prad_tcga_pub_sequenced"
)
save(cna_tcga, file="data-raw/cna_tcga.RData")

# Mutations
mut_tcga <- curatedPCaData:::generate_cbioportal(
  genes = sort(unique(curatedPCaData:::curatedPCaData_genes$hgnc_symbol)),
  geneticProfiles = "prad_tcga_pub_mutations",
  caseList="prad_tcga_pub_sequenced"
)
# Save NA values as truly NA instead of "NaN" even if other instances exist on column
mut_tcga[which(mut_tcga=="NaN")] <- NA
save(mut_tcga, file="data-raw/mut_tcga.RData")

# Create MAE object
mae_tcga <- curatedPCaData:::create_mae(study_name = "TCGA")
usethis::use_data(mae_tcga, overwrite = TRUE)

## - end TCGA -

## - True et al. -
# GEX:
#	GPL3834	FHCRC Human Prostate PEDB cDNA Array v4
#	GPL3836	FHCRC Human Prostate PEDB cDNA Array v3
gex_true <- curatedPCaData:::generate_gex_geo(
	geo_code = "GSE5132",
	pckg = "oligo",
	cleanup = FALSE
)
save(gex_true, file = "data-raw/gex_true.RData")

# Create MAE object
mae_true <- curatedPCaData:::create_mae(study_name = "true")
usethis::use_data(mae_true, internal = FALSE, overwrite = TRUE)

## - end True et al. -


## - Wallace et al. -
# GEX: [HG-U133A_2] Affymetrix Human Genome U133A 2.0 Array
gex.rma_wallace <- curatedPCaData:::generate_gex_geo(
	geo_code = "GSE6956",
	pckg = "oligo",
	cleanup = FALSE
)
save(gex.rma_wallace, file = "data-raw/gex.rma_wallace.RData")

# Create and save MAE object
mae_wallace <- curatedPCaData:::create_mae(study_name = "wallace")
usethis::use_data(mae_wallace, internal = FALSE, overwrite = TRUE)

## - end Wallace et al. -


# - Wang et al. -
# GEX: [HG-U133A] Affymetrix Human Genome U133A Array
gex.rma_wang <- curatedPCaData:::generate_gex_geo(
	geo_code = "GSE8218",
	pckg = "oligo"
)
save(gex.rma_wang, file="data-raw/gex.rma_wang.RData")
#CNA
cna_wang <- curatedPCaData:::generate_cna_geo(
  geo_code = "GSE8218"
)

save(cna_wang, file="data-raw/cna_wang.RData")

# Create MAE object
mae_wang <- curatedPCaData:::create_mae(study_name = "wang")
usethis::use_data(mae_wang, overwrite = TRUE)

# - end Wang et al. -








