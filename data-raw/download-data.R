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
# Friedrich 2020 BASED ON GPL DATA
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

##
#
# Chandran et al.
#
##

# Create and save GEX of Chandran et al., Yu et al.
gex_chandran <- curatedPCaData::generate_gex_geo("GSE6919")
save(gex_chandran, file="data-raw/gex_chandran.RData")

# Create and save MAE object
mae_chandran <- curatedPCaData:::create_mae(study_name = "chandran")
usethis::use_data(mae_chandran, internal = FALSE, overwrite = TRUE)

######################################################################
#Wallace et a. BASED ON RAW AFFY DATA
######################################################################

library(hgu133a2cdf, lib = '~/Rloc')
library(hgu133a.db, lib = '~/Rloc')
collapseFUN = function(z) { apply(z, MARGIN=2, FUN=median) }

unmatched_healty_tissue = c('GSM160418', 'GSM160419', 'GSM160420', 'GSM160421', 'GSM160422', 'GSM160430')



supfiles <- GEOquery::getGEOSuppFiles('GSE6956')

utils::untar(tarfile=rownames(supfiles))
	# Make sure to function in a working directory where the are no other tarballs present
supfiles2 <- list.files()
supfiles2 <- supfiles2[grep(".gz", supfiles2)]

Wallace <- affy::ReadAffy()
colnames(affy::exprs(Wallace)) <- gsub(".gz|.CEL", "", colnames(Wallace))
GEX_Wallace <- affy::rma(Wallace)
# Removing .CEL and packaging names from the GEO-compatible sample names
colnames(GEX_Wallace) <- gsub(".CEL.gz", "", colnames(affy::exprs(GEX_Wallace)))
keys <- AnnotationDbi::mappedkeys(hgu133a.db::hgu133aGENENAME)
nam <- names(as.character(hgu133a.db::hgu133aALIAS2PROBE)[match(rownames(GEX_Wallace), as.character(hgu133a.db::hgu133aALIAS2PROBE))])
nam[is.na(nam)] <- "NA"
# Collapse probes
gex_wallace <- do.call("rbind", by(as.matrix(affy::exprs(GEX_Wallace)), INDICES=nam, FUN=collapseFUN))
# Remove downloaded files
if(cleanup){
  # First GEO download
  file.remove(rownames(supfiles))
  # Tarballs
  file.remove(supfiles2)
}
# Return numeric matrix

gex_wallace = gex_wallace[, !is.element(colnames(gex_wallace), unmatched_healty_tissue)]

save(gex_wallace,  file = "./gex_wallace.RData")


#Download OSF data
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
mae_ren <- create_mae(study_name = "ren")
usethis::use_data(mae_ren, overwrite = TRUE)

#Kim et al data -----
gex_kim <- generate_gex_geo(
  geo_code = "GSE119616"
)
save(gex_kim, file="data-raw/gex_kim.RData")

# CNA
cna_kim <- generate_cna_geo(
  geo_code = "GSE119616"
)
save(cna_kim, file="data-raw/cna_kim.RData")

# Create MAE object
mae_kim <- create_mae(study_name = "kim")
usethis::use_data(mae_kim, overwrite = TRUE)

#Abida et al data -----
# GEX Capture assay (FPKM)
gex_capture_abida <- generate_cbioportal(
  genes = sort(unique(curatedPCaData_genes$hgnc_symbol)), # All unique gene symbols
  geneticProfiles = "prad_su2c_2019_mrna_seq_fpkm_capture_all_sample_Zscores", # Omics profile
  caseList = "prad_su2c_2019_all" # Case list
)
save(gex_capture_abida, file="data-raw/gex_capture_abida.RData")

# GEX PolyA (FPKM)
gex_polyA_abida <- generate_cbioportal(
  genes = sort(unique(curatedPCaData_genes$hgnc_symbol)), # All unique gene symbols
  geneticProfiles = "prad_su2c_2019_mrna_seq_fpkm_polya_all_sample_Zscores", # Omics profile
  caseList = "prad_su2c_2019_all" # Case list
)
save(gex_polyA_abida, file="data-raw/gex_polyA_abida.RData")

# CNA
cna_abida <- generate_cbioportal(
  genes = sort(unique(curatedPCaData_genes$hgnc_symbol)),
  geneticProfiles="prad_su2c_2019_gistic", 
  caseList="prad_su2c_2019_sequenced"
)
save(cna_abida, file="data-raw/cna_abida.RData")

# Create MAE object
mae_abida <- create_mae(study_name = "abida")
usethis::use_data(mae_abida, overwrite = TRUE)

# Wang et al.
#GEX
gex_wang <- curatedPCaData:::generate_gex_geo(
  geo_code = "GSE8218"
)
save(gex_wang, file="data-raw/gex_wang.RData")
#CNA
cna_wang <- curatedPCaData:::generate_cna_geo(
  geo_code = "GSE8218"
)
save(cna_wang, file="data-raw/cna_wang.RData")

# Barwick et al.
gex_barwick <- generate_gex_geo("GSE18655")
save(gex_barwick, file="data-raw/gex_barwick.RData")


# Kunderfranco

g_kunder <- getGEO("GSE14206", GSEMatrix =TRUE, getGPL=TRUE)
kunder_labels = Biobase::fData(g_kunder[[1]])
kunder_ex <- exprs(g_kunder[[1]])
rownames(kex) = kunder_labels$GENE_SYMBOL
kunder_ex = kunder_ex[-which(rownames(kunder_ex) == ''), ]

gex_kunderfranco = aggregate(kunder_ex, by = list(rownames(ex)), mean)

rownames(gex_kunderfranco) = gex_kunderfranco[, 1]
gex_kunderfranco = gex_kunderfranco[, -1]

save(gex_kunderfranco, file = "data-raw/gex_kunderfranco.RData")



