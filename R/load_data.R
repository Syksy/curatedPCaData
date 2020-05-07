# Paper to read
# http://tau.amegroups.com/article/view/24493/23261

# Load data

############################################################################## TCGA ############# 
# RTCGA: The Cancer Genome Atlas Data Integration # https://rdrr.io/bioc/RTCGA/man/downloadTCGA.html
# TCGAretriever helps accessing and downloading TCGA data hosted on 'cBioPortal' via its Web Interface
# This one looks better
# TCGAbiolinks package # https://bioc.ism.ac.jp/packages/3.2/bioc/vignettes/TCGAbiolinks/inst/doc/tcgaBiolinks.html

# install packages
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.10")
# library(BiocManager)
# BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)

# Load data
TCGAbiolinks:::getProjectSummary("TCGA-PRAD")
# $data_categories
# case_count file_count               data_category
# 1        498       2755     Transcriptome Profiling
# 2        498       4032 Simple Nucleotide Variation
# 3        500       2182                 Biospecimen
# 4        500        537                    Clinical
# 5        498        553             DNA Methylation
# 6        498       3059       Copy Number Variation
# 7        498       2178            Sequencing Reads
query <- GDCquery(project = "TCGA-PRAD", 
                  data.category = "Clinical",
                  data.type = "Clinical Supplement", 
                  data.format = "BCR Biotab"
                  )
clinical_tab_all <- GDCdownload(query)
clinical_tab_all <- GDCprepare(query)
clinical_tcga <- Reduce(function(...) merge(..., by="bcr_patient_barcode", all=TRUE), clinical_tab_all) 
query <- GDCquery(project = "TCGA-PRAD", 
                  data.category = "Clinical",
                  data.type = "Clinical Supplement", 
                  data.format = "BCR XML"
)
clinical <- GDCprepare_clinic(query, clinical.info = "patient")
clin <- merge.data.frame(clinical_tcga, clinical, 
                 by.x = "bcr_patient_barcode", by.y = "bcr_patient_barcode",
                 all.x = TRUE, all.y = TRUE)
# clinical_tcga$American Joint Committee on Cancer Tumor Stage Code
# 
# clinical_tab_al <- GDCprepare(query)
# clinXML <- bind_rows(clinical_tab_al)
# 
# query <- GDCquery(project = "TCGA-PRAD", 
#                   data.category = "Biospecimen",
#                   data.type = "Biospecimen Supplement", 
#                   data.format = "BCR Biotab"
# )
# clinical_tab_all <- GDCdownload(query)
# clinical_tab_all <- GDCprepare(query)
# clinical_tab_all$ssf_tumor_samples_prad
# 
# b <- merge(clinical_tab_all$ssf_normal_controls_prad, clinical_tab_all$biospecimen_diagnostic_slides_prad)
# colnames(b)
# c <- merge(b, clinical_tab_all$biospecimen_slide_prad)
# d <- merge.data.frame(c, clinical_tab_all$ssf_tumor_samples_prad, 
#                       by.x = "bcr_patient_barcode", by.y = "bcr_patient_barcode",
#                       all.x = TRUE, all.y = TRUE)
# e <- merge(d, clinical_tab_all$biospecimen_analyte_prad)
# f <- merge(e, clinical_tab_all$biospecimen_sample_prad)

# a <- merge(clinical_tcga, f)
# colnames(a)


query1 <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Clinical", 
                  file.type = "xml", 
                  barcode = c("TCGA-RU-A8FL","TCGA-AA-3972"))
GDCdownload(query1)


clinical_tcga <- Reduce(function(...) merge(..., by="bcr_patient_barcode", all=TRUE), clinical_tab_all) 
colnames(clinical_tcga)
#clinical_tcga$last_contact_days_to.x == clinical_tcga$last_contact_days_to.y

# Cleaning
rm(query, clinical_tab_all)
############################################################################## From cBioPortal 

############################################################################## SUN ############# 
# Use package cgdsr
# install.packages('cgdsr')
library(cgdsr)
# Create CGDS object
mycgds = CGDS("https://www.cbioportal.org/")
test(mycgds)

# Get list of cancer studies at server
# a <- getCancerStudies(mycgds) # We want prad_eururol_2017

# Get available case lists (collection of samples) for a given cancer study
# mycaselist = getCaseLists(mycgds,"prad_eururol_2017") # In here I want prad_eururol_2017_all

# Get clinical data for the case list
clinical_SUN = getClinicalData(mycgds,"prad_eururol_2017_all")


############################################################################## TCGA_333 ############# 
# Get clinical data for the case list
clinical_TCGA_333 = getClinicalData(mycgds,"prad_tcga_pub_all")


clinical_TCGA_333$bcr_patient_barcode <- substr(rownames(clinical_TCGA_333), start = 1, stop = 12)

a <- merge.data.frame(clinical_TCGA_333, clinical[, c("bcr_patient_barcode", "stage_event_tnm_categories")], 
                      by.x = "bcr_patient_barcode", by.y = "bcr_patient_barcode",
                      all.x = TRUE, all.y = TRUE)
# clinical_TCGA_333b <- a %>% 
#   select(c("bcr_patient_barcode", "PREOPERATIVE_PSA", "REVIEWED_GLEASON", "REVIEWED_GLEASON_CATEGORY",
#            "REVIEWED_GLEASON_SUM", "AGE", "clinical_T"))
colnames(a)
table(a$stage_event_tnm_categories)
############################################################################## MSKCC ############# 
clinical_MSKCC = getClinicalData(mycgds,"prad_cdk12_mskcc_2020_all")


############################################################################## CPC-GENE ############# 
clinical_CPC_GENE = getClinicalData(mycgds,"prad_cpcg_2017_all")


############################################################################## From ICGC 

############################################################################## ICGC-FR #############
icgc <- list()
icgc[["PRAD-FR"]] <- 
  c(
    "https://dcc.icgc.org/api/v1/download?fn=/current/Projects/PRAD-FR/donor.tsv.gz")

install.packages("remotes")
remotes::install_github("Syksy/curatedTools")
clinical_ICGC_FR <- lapply(icgc[["PRAD-FR"]], FUN=curatedTools:::.icgcDownload)

icgc[["PRAD-UK"]] <- 
  c(
    "https://dcc.icgc.org/api/v1/download?fn=/current/Projects/PRAD-UK/donor.PRAD-UK.tsv.gz")

