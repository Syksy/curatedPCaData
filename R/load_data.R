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
query <- GDCquery(project = "TCGA-PRAD", 
                  data.category = "Clinical",
                  data.type = "Clinical Supplement", 
                  data.format = "BCR Biotab"
                  )
clinical_tab_all <- GDCprepare(query)


clinical_tcga <- Reduce(function(...) merge(..., by="bcr_patient_barcode", all=TRUE), clinical_tab_all) #%>% 
# select()

colnames(clinical_tcga)
clinical_tcga$last_contact_days_to.x == clinical_tcga$last_contact_days_to.y
table(clinical_tcga$therapy_regimen.y)
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
a <- merge.data.frame(clinical_TCGA_333, clinical_tcga, 
                      by.x = "bcr_patient_barcode", by.y = "bcr_patient_barcode",
                      all.x = TRUE, all.y = FALSE)
# clinical_TCGA_333b <- a %>% 
#   select(c("bcr_patient_barcode", "PREOPERATIVE_PSA", "REVIEWED_GLEASON", "REVIEWED_GLEASON_CATEGORY",
#            "REVIEWED_GLEASON_SUM", "AGE", "clinical_T"))
colnames(a)

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
ICGC.PRAD.FR <- lapply(icgc[["PRAD-FR"]], FUN=curatedTools:::.icgcDownload)
clinical_ICGC_FR <- lapply(icgc[["PRAD-FR"]], FUN=curatedTools:::.icgcDownload)

icgc[["PRAD-UK"]] <- 
  c(
    "https://dcc.icgc.org/api/v1/download?fn=/current/Projects/PRAD-UK/donor.PRAD-UK.tsv.gz")
