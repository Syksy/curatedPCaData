# Import library

library(tidyverse)

# Paper to read
# http://tau.amegroups.com/article/view/24493/23261

# Load data

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
clinical_TCGA_333 = getClinicalData(mycgds,"prad_tcga_pub_all")
# Rename the patients id for later merge
clinical_TCGA_333$bcr_patient_barcode <- substr(rownames(clinical_TCGA_333), start = 1, stop = 12) %>% 
  str_replace_all("\\.", "-")

# Get clinical data for the case list from biolinks to have clinical T stage

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
                  data.format = "BCR XML"
)
clinical <- GDCprepare_clinic(query, clinical.info = "patient")
# Merge and keep only the 333 patients
clinical_TCGA_333 <- merge.data.frame(
  clinical_TCGA_333, 
  clinical[, c("bcr_patient_barcode", "stage_event_tnm_categories")], 
                      by.x = "bcr_patient_barcode", by.y = "bcr_patient_barcode",
                      all.x = TRUE, all.y = FALSE)


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

# Cleaning
rm(query, clinical, icgc, mycgds)
