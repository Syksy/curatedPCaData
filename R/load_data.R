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
############################################################################## For cBioPortal 

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

############################################################################## MSKCC ############# 
# Get clinical data for the case list
clinical_MSKCC = getClinicalData(mycgds,"prad_cdk12_mskcc_2020_all")

############################################################################## CPC-GENE ############# 
# Get clinical data for the case list
clinical_CPC_GENE = getClinicalData(mycgds,"prad_cpcg_2017_all")

