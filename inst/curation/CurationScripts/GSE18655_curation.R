rm(list=ls())
source("functions.R")

uncurated <- read.delim("../DFCI_version/uncurated_copy/GSE18655_full_pdata.csv",as.is=TRUE,row.names=1,sep=",",header=T)

##initial creation of curated dataframe
curated <- initialCuratedDF(rownames(uncurated),template.filename="template_prad.csv")

##--------------------
##start the mappings
##--------------------

# two additional categories that might be useful, but not sure what to classify them as
# pathologic stage category
# grade (1-4)

##study name
curated$study_name <- rep("Barwick, et al.", dim(curated)[1])

##alternative sample name and patient id
curated$patient_id <- uncurated$geo_accession
curated$alt_sample_name <- uncurated$title
  
#age at diagnosis
tmp <- uncurated$characteristics_ch1
tmp <- gsub("age: ", "",tmp)
curated$age_at_initial_diagnosis <- as.numeric(tmp)

#tumor sample source
tmp <- uncurated$source_name_ch1
tmp <- gsub("Radical ", "",tmp)
curated$tissue_source <- tmp

#PSA
tmp <- uncurated$characteristics_ch1.1
tmp <- gsub("psa: ", "",tmp)
curated$psa <- as.numeric(tmp)

#gleason grade
tmp <- uncurated$characteristics_ch1.2
tmp <- gsub("gleason score: ", "",tmp)
curated$gleason_grade <- as.numeric(tmp)

##disease free survival status
tmp <- uncurated$characteristics_ch1.3
tmp <- gsub("recurrence: ", "",tmp)
curated$disease_specific_recurrence_status <- ifelse(grepl("No Rec",tmp,ignore.case=TRUE),0,1)

##disease free survival months
tmp <- uncurated$characteristics_ch1.4
tmp <- gsub("follow-up \\(months\\): ", "",tmp)
tmp <- as.numeric(tmp) * 30.5
curated$days_to_disease_specific_recurrence <- tmp 

##ERG fusion status
tmp <- uncurated$characteristics_ch1.5
tmp <- gsub("tmprss2: ERG Fusion: ", "",tmp)
curated$ERG_fusion_CNA <- ifelse(grepl("No Fusion",tmp,ignore.case=TRUE),0,1) # called from PCR

##Positive margins
tmp <- uncurated$characteristics_ch1.6
tmp <- gsub("positive surgical margin: ", "",tmp)
curated$tumor_margins_positive <- ifelse(grepl("No Positive Margin",tmp,ignore.case=TRUE),0,1)


write.table(curated, row.names = FALSE, file="../DFCI_version/curated/GSE18655_curated_pdata.txt",sep="\t")

