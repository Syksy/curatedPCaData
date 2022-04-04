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
## - Abida et al.
## - Baca et al.
## - Barbieri et al.
## - Barwick et al. 
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
## - Weiner et al.
## (to be updated)

## - Abida et al. -
# GEX PolyA (FPKM) z-score normalized relative to paired normal
# gex.relz_abida <- curatedPCaData:::generate_cbioportal(
#   genes = sort(unique(curatedPCaData:::curatedPCaData_genes$hgnc_symbol)), # All unique gene symbols
#   geneticProfiles = "prad_su2c_2019_mrna_seq_fpkm_polya_all_sample_Zscores", # Omics profile
#   caseList = "prad_su2c_2019_all" # Case list
# )
gex.relz_abida <- curatedPCaData:::generate_cbioportaldata("prad_su2c_2019","gex")
save(gex.relz_abida, file="data-raw/gex.relz_abida.RData")

# CNA (discretized GISTIC)
# cna.gistic_abida <- curatedPCaData:::generate_cbioportal(
#   genes = sort(unique(curatedPCaData:::curatedPCaData_genes$hgnc_symbol)),
#   geneticProfiles="prad_su2c_2019_gistic", 
#   caseList="prad_su2c_2019_sequenced"
# )
cna.gistic_abida <- curatedPCaData:::generate_cbioportaldata("prad_su2c_2019","cna")
save(cna.gistic_abida, file="data-raw/cna.gistic_abida.RData")

# Mutations
# mut_abida <- curatedPCaData:::generate_cbioportal(
#   genes = sort(unique(curatedPCaData:::curatedPCaData_genes$hgnc_symbol)),
#   geneticProfiles="prad_su2c_2019_mutations", 
#   caseList="prad_su2c_2019_sequenced"
# )
# mut_abida[which(mut_abida=="NaN")] <- NA
abida_mut <- curatedPCaData:::generate_cbioportaldata("prad_su2c_2019","mut")
save(abida_mut, file="data-raw/mut_abida.RData")
# To check: Fusions separately?

# Create MAE object
mae_abida <- curatedPCaData:::create_mae(study_name = "abida")
usethis::use_data(mae_abida, overwrite = TRUE)

## - end Abida et al. -


## - Baca et al. -
# CNA
# cna.gistic_baca <- curatedPCaData:::generate_cbioportal(
#   genes = sort(unique(curatedPCaData:::curatedPCaData_genes$hgnc_symbol)),
#   geneticProfiles="prad_broad_2013_cna", 
#   caseList="prad_broad_2013_sequenced"
# )
cna.gistic_baca <-curatedPCaData:::generate_cbioportaldata("prad_broad_2013","cna")
save(cna.gistic_baca, file="data-raw/cna.gistic_baca.RData")

# Mutations
baca_mut <- curatedPCaData:::generate_cbioportaldata("prad_broad_2013","mut")
save(baca_mut, file="data-raw/mut_baca.RData")

# Create MAE object
mae_baca <- curatedPCaData:::create_mae(study_name = "baca")
usethis::use_data(mae_baca, overwrite = TRUE)

## - end Baca et al. -


## - Barbieri et al. -
# GEX expression z-score normalized relative to paired normal
# gex.relz_barbieri <- curatedPCaData:::generate_cbioportal(
#   genes = sort(unique(curatedPCaData:::curatedPCaData_genes$hgnc_symbol)), # All unique gene symbols
#   geneticProfiles = "prad_broad_mrna", # Omics profile
#   caseList = "prad_broad_sequenced" # Case list
# )
gex.relz_barbieri <- curatedPCaData:::generate_cbioportaldata("prad_broad","gex")
save(gex.relz_barbieri, file="data-raw/gex.relz_barbieri.RData")

# CNA
# cna.gistic_barbieri <- curatedPCaData:::generate_cbioportal(
#   genes = sort(unique(curatedPCaData:::curatedPCaData_genes$hgnc_symbol)),
#   geneticProfiles="prad_broad_cna", 
#   caseList="prad_broad_sequenced"
# )
cna.gistic_barbieri <- curatedPCaData:::generate_cbioportaldata("prad_broad","cna")
save(cna.gistic_barbieri, file="data-raw/cna.gistic_barbieri.RData")

# MUT
# mut_barbieri <- curatedPCaData:::generate_cbioportal(
#   genes = sort(unique(curatedPCaData:::curatedPCaData_genes$hgnc_symbol)),
#   geneticProfiles="prad_broad_mutations", 
#   caseList="prad_broad_sequenced"
# )
# mut_barbieri[which(mut_barbieri=="NaN")] <- NA

# barbieri_mut<-curatedPCaData:::generate_cbioportaldata_mut(
#   caselist = "prad_broad"
#   )
barbieri_mut <- curatedPCaData:::generate_cbioportaldata("prad_broad","mut")
save(barbieri_mut, file="data-raw/mut_barbieri.RData")

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
gex.rma_chandran <- curatedPCaData::generate_gex_geo(
	geo_code = "GSE6919",
	pckg = "oligo",
	filter_regex = "_RAW"
)
save(gex.rma_chandran, file="data-raw/gex.rma_chandran.RData")

# Create and save MAE object
mae_chandran <- curatedPCaData:::create_mae(study_name = "chandran")
usethis::use_data(mae_chandran, internal = FALSE, overwrite = TRUE)

# - end Chandran et al. -


## - Friedrich et al. (2020) -
# GEX: GPL26898	Agilent-058029 Custom human expression microarray (Probe Name version)
gex.logq_friedrich <- curatedPCaData:::generate_gex_geo(
	geo_code = "GSE134051",
	pckg = "limma",
	filter_regex = "_RAW",
	cleanup = FALSE
)
save(gex.logq_friedrich, file="data-raw/gex.logq_friedrich.RData")

# Create and save MAE object
mae_friedrich = curatedPCaData:::create_mae(study_name = 'Friedrich')
usethis::use_data(mae_friedrich, internal = FALSE, overwrite = TRUE)

## - end Friedrich et al. -


## - Hieronymus et al. -
# CNA
cna.logr_hieronymus <- curatedPCaData:::generate_cna_geo(
  geo_code = "GSE54691"
)
save(cna.logr_hieronymus, file="data-raw/cna.logr_hieronymus.RData")

# Create MAE object
mae_hieronymus <- curatedPCaData:::create_mae(study_name = "hieronymus")
usethis::use_data(mae_hieronymus, internal = FALSE, overwrite = TRUE)

## - end Hieronymus et al. -


## - ICGC datasets -

# PRAD-CA
# Affymetrix Human Gene 1.0 ST
gex.rma_icgcca <- curatedPCaData:::generate_icgc("PRAD_CA", "gex")
save(gex.rma_icgcca, file="data-raw/gex.rma_icgcca.RData")

# TODO: Copy-number alterations, mutations (available in ICGC)

# Create MAE object
mae_icgcca <- curatedPCaData:::create_mae(study_name = "icgcca")
usethis::use_data(mae_icgcca, internal = FALSE, overwrite = TRUE)



# PRAD-FR
gex_icgcfr <- curatedPCaData::generate_icgc("PRAD_FR", "gex")
save(gex_icgcfr, file="data-raw/gex_icgcfr.RData")

# TODO: At this point the data contains raw read counts, and is not yet usable as a 2-dim gex matrix

# TODO: MAE

# PRAD-UK

# Only contains CNA data

## - end ICGC datatasets -



## - IGC -
# GEX: [HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array
gex.rma_igc <- curatedPCaData:::generate_gex_geo(
	geo_code = "GSE2109",
	cleanup = FALSE,
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
	geo_code = "GSE119616",
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
# Global LOESS normalized log ratios between the two colour arrays combined by mean with flipped dye swap
gex.logr_kunderfranco <- curatedPCaData:::generate_gex_geo(
	geo_code = "GSE14206",
	pckg = "limma",
	cleanup = FALSE	
)
save(gex.logr_kunderfranco, file = "data-raw/gex.logr_kunderfranco.RData")

# Create MAE object
mae_kunderfranco <- curatedPCaData:::create_mae(study_name = "kunderfranco")
usethis::use_data(mae_kunderfranco, internal = FALSE, overwrite = TRUE)

## - end Kunderfranco et al. -



## - Ren et al. -
# GEX
# gex.relz_ren <- curatedPCaData:::generate_cbioportal(
#   genes = sort(unique(curatedPCaData:::curatedPCaData_genes$hgnc_symbol)), # All unique gene symbols
#   geneticProfiles = "prad_eururol_2017_rna_seq_mrna", # Omics profile
#   caseList = "prad_eururol_2017_sequenced" # Case list
# )
gex.relz_ren <- curatedPCaData:::generate_cbioportaldata("prad_eururol_2017","gex")
save(gex.relz_ren, file="data-raw/gex.relz_ren.RData")

# CNA
# cna.gistic_ren <- curatedPCaData:::generate_cbioportal(
#   genes = sort(unique(curatedPCaData:::curatedPCaData_genes$hgnc_symbol)),
#   geneticProfiles="prad_eururol_2017_cna", 
#   caseList="prad_eururol_2017_sequenced"
# )
cna.gistic_ren <- curatedPCaData:::generate_cbioportaldata("prad_eururol_2017","cna")
save(cna.gistic_ren, file="data-raw/cna.gistic_ren.RData")

# Mutations
# mut_ren <- curatedPCaData:::generate_cbioportal(
#   genes = sort(unique(curatedPCaData:::curatedPCaData_genes$hgnc_symbol)), # All unique gene symbols
#   geneticProfiles = "prad_eururol_2017_mutations", # Omics profile
#   caseList = "prad_eururol_2017_sequenced" # Case list
# )
# mut_ren[which(mut_ren=="NaN")] <- NA
# ren_mut<-curatedPCaData:::generate_cbioportaldata_mut(
#   caselist = "prad_eururol_2017"
#   )
ren_mut <- curatedPCaData:::generate_cbioportaldata("prad_eururol_2017","mut")
save(ren_mut, file="data-raw/mut_ren.RData")

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
save(gex.rma_taylor, file="data-raw/gex.rma_taylor.RData")

# CNA
cna.logr_taylor <- curatedPCaData:::generate_cna_geo(
	geo_code = "GSE21035"
)
save(cna.logr_taylor, file="data-raw/cna.logr_taylor.RData")

cna.gistic_taylor <- curatedPCaData:::generate_cbioportaldata(
  "prad_mskcc","cna"
)
save(cna.gistic_taylor, file="data-raw/cna.gistic_taylor.RData")

# Mutations - notice this is downloaded from cBioPortal rather than processed from GEO
# mut_taylor <- curatedPCaData:::generate_cbioportal(
# 	genes = sort(unique(curatedPCaData:::curatedPCaData_genes$hgnc_symbol)),
# 	geneticProfiles = "prad_mskcc_mutations",
# 	caseList="prad_mskcc_sequenced"
# ) 
# # Save NA values as truly NA instead of "NaN" even if other instances exist on column
# mut_taylor[which(mut_taylor=="NaN")] <- NA
# # Grep down to using only patient samples, omitting cell lines etc
# mut_taylor <- mut_taylor[,grep("PCA", colnames(mut_taylor))]
taylor_mut<- curatedPCaData:::generate_cbioportaldata("prad_mskcc","mut")

save(taylor_mut, file="data-raw/mut_taylor.RData")

# Create MAE object
mae_taylor <- curatedPCaData:::create_mae(study_name = "Taylor")
usethis::use_data(mae_taylor, overwrite = TRUE)

## - end Taylor et al. (MSKCC) -


## - TCGA - 
# GEX
gex.fpkm_tcga <- curatedPCaData:::generate_xenabrowser(
	id = "TCGA-PRAD",
	type = "gex",
	digits = 4,
	truncate = 0 # '.01A' -> '.01' suffix
)
save(gex.fpkm_tcga, file="data-raw/gex.fpkm_tcga.RData")

# CNA (GISTIC)
cna.gistic_tcga <- curatedPCaData:::generate_xenabrowser(
	id = "TCGA-PRAD",
	type = "cna",
	truncate = 0 # '.01A' -> '.01' suffix
)
save(cna.gistic_tcga, file="data-raw/cna.gistic_tcga.RData")

# Mutations
tcga_mut <- curatedPCaData:::generate_xenabrowser(
  id = "TCGA-PRAD",
  type = "mut",
  truncate = 0 # '.01A' -> '.01' suffix
)
# Save NA values as truly NA instead of "NaN" even if other instances exist on column
#mut_tcga[which(mut_tcga=="NaN")] <- NA
save(tcga_mut, file="data-raw/mut_tcga.RData")

# Create MAE object
mae_tcga <- curatedPCaData:::create_mae(study_name = "TCGA")
usethis::use_data(mae_tcga, overwrite = TRUE)

## - end TCGA -

## - True et al. -
# GEX:
#	GPL3834	FHCRC Human Prostate PEDB cDNA Array v4
#	GPL3836	FHCRC Human Prostate PEDB cDNA Array v3 (-> single sample only! 11th, GSM115769)
gex.logr_true <- curatedPCaData:::generate_gex_geo(
	geo_code = "GSE5132",
	pckg = "limma",
	cleanup = FALSE,
	filter_regex = "_RAW"
)
save(gex.logr_true, file = "data-raw/gex.logr_true.RData")

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

##CNA
#cna_wang <- curatedPCaData:::generate_cna_geo(
#  geo_code = "GSE8218"
#)
#save(cna_wang, file="data-raw/cna_wang.RData")

# Create MAE object
mae_wang <- curatedPCaData:::create_mae(study_name = "wang")
usethis::use_data(mae_wang, overwrite = TRUE)

# - end Wang et al. -


# - Weiner et al. -
# GEX: GPL5175	[HuEx-1_0-st] Affymetrix Human Exon 1.0 ST Array [transcript (gene) version]
gex.rma_weiner <- curatedPCaData:::generate_gex_geo(
	geo_code = "GSE157548",
	pckg = "oligo"
)
save(gex.rma_weiner, file="data-raw/gex.rma_weiner.RData")

# Create MAE object
mae_weiner <- curatedPCaData:::create_mae(study_name = "weiner")
usethis::use_data(mae_weiner, overwrite = TRUE)

# - end Weiner et al. -





