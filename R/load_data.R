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
# rm(clinical_tab_all)
colnames(clinical_tcga)
clinical_tcga$last_contact_days_to.x == clinical_tcga$last_contact_days_to.y
table(clinical_tcga$therapy_regimen.y)
