rm(list=ls())
source("functions.R")

uncurated <- read.delim("../DFCI_version/uncurated_copy/GSE16560_full_pdata.csv",as.is=TRUE,row.names=1,sep=",")

##initial creation of curated dataframe
##curated <- data.frame(sample_name=rownames(uncurated))#,row.names=rownames(uncurated))
##initial creation of curated dataframe
curated <- initialCuratedDF(rownames(uncurated),template.filename="template_prad.csv")

##--------------------
##start the mappings
##--------------------

##study name
curated$study_name <- rep("Sboner, et al.", dim(curated)[1])

##alternative sample name and patient id
tmp <- uncurated$title
curated$alt_sample_name <- tmp
curated$patient_id <- tmp


#tumor sample site
tmp <- uncurated$source_name_ch1
tmp <- gsub(" prostate cancer", "",tmp)
curated$sample_type <- tmp

#gleason grade
tmp <- uncurated$characteristics_ch1
tmp <- gsub("gleason: ", "",tmp)
tmp <- as.numeric(tmp)
curated$gleason_grade <- tmp

#gleason major
tmp <- uncurated$characteristics_ch1.1
tmp <- gsub("major.gleason: ", "",tmp)
tmp <- as.numeric(tmp)
curated$gleason_major <- tmp

#gleason minor
tmp <- uncurated$characteristics_ch1.2
tmp <- gsub("minor.gleason: ", "",tmp)
tmp <- as.numeric(tmp)
curated$gleason_minor <- tmp

#gleason group
tmp <- uncurated$characteristics_ch1
tmp <- gsub("gleason: ", "",tmp)
tmp <- gsub("6", "<=6", tmp)
tmp <- gsub("8", ">=8", tmp)
tmp <- gsub("9", ">=8", tmp)
tmp <- gsub("10", ">=8", tmp)
tmp[curated$gleason_major == 3 & curated$gleason_minor == 4] <- "3+4"
tmp[curated$gleason_major == 4 & curated$gleason_minor == 3] <- "4+3"
curated$grade_group

#batch
tmp <- uncurated$characteristics_ch1.3
tmp <- gsub("batch: ", "",tmp)
tmp <- as.numeric(tmp)
curated$batch <- tmp

#tumor purity pathology 
tmp <- uncurated$characteristics_ch1.5
tmp <- gsub("cancer.percent: ", "",tmp)
tmp <- as.numeric(tmp)
curated$tumor_purity_pathology <- tmp

#age at diagnosis
tmp <- uncurated$characteristics_ch1.6
tmp <- gsub("age: ", "",tmp)
tmp <- as.numeric(tmp)
curated$age_at_initial_diagnosis <- tmp 

#year at diagnosis
tmp <- uncurated$characteristics_ch1.7
tmp <- gsub("diag.yr: ", "",tmp)
tmp <- as.numeric(tmp)
curated$year_diagnosis <- tmp 

##overall survival status
tmp <- uncurated$characteristics_ch1.8
curated$overall_survival_status <- ifelse(grepl("Dead",tmp,ignore.case=TRUE),1,0)
curated$overall_survival_status[is.na(tmp)] <- NA

##overall survival months
tmp <- uncurated$characteristics_ch1.9
tmp <- gsub("fup.month: ", "",tmp)
tmp <- as.numeric(tmp)
tmp <- tmp * 30.5
curated$days_to_overall_survival <- tmp 

##ERG fusion status
tmp <- uncurated$characteristics_ch1.10
tmp <- gsub("fusion: ", "",tmp)
tmp <- as.numeric(tmp)
curated$ERG_fusion_IHC <- tmp 

write.table(curated, row.names = FALSE, file="../DFCI_version/curated/GSE16560_curated_pdata.txt",sep="\t")
