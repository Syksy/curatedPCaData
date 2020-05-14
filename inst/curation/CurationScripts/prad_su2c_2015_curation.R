rm(list=ls())
source("functions.R")

uncurated <- read.delim("../DFCI_version/uncurated_copy/prad_su2c_2015/data_clinical_patient.txt",as.is=TRUE,row.names=1,sep="\t", skip=4)
colnames(uncurated) <- as.character(read.delim("../DFCI_version/uncurated_copy/prad_su2c_2015/data_clinical_patient.txt",as.is=TRUE,header=T,row.names=1,sep="\t",skip=3)[1,])

##initial creation of curated dataframe
##curated <- data.frame(sample_name=rownames(uncurated))#,row.names=rownames(uncurated))
##initial creation of curated dataframe
curated <- initialCuratedDF(rownames(uncurated),template.filename="template_prad.csv")

##--------------------
##start the mappings
##--------------------

####################################
# some treatment information might be included

##study name
curated$study_name <- rep("Robinson, et al.", dim(curated)[1])

##sample type is all mets in this study
curated$sample_type <- rep("metastatic", dim(curated)[1])

##patient id
curated$patient_id <- row.names(uncurated)
curated$sample_name <- row.names(uncurated)

#age at diagnosis
curated$age_at_initial_diagnosis <- uncurated$AGE 

#hormonal treatment status (abi and enza status reported)
tmp <- uncurated$ABI_ENZA_EXPOSURE_STATUS
tmp <- gsub("Yes", "1",tmp)
tmp <- gsub("No", "0",tmp)
tmp <- gsub("", "NA",tmp)
curated$therapy_hormonal_initial <- as.numeric(tmp)

#hormonal treatment status (abi and enza status reported)
tmp <- uncurated$PRIOR_TREATMENT
tmp <- ifelse(grepl("Had prior Taxane treatment",tmp,ignore.case=TRUE),"taxane","")
curated$other_treatment <- tmp

#metastatic tumor site
tmp <- uncurated$TUMOR_SITE
tmp <- gsub("Bone", "bone",tmp)
tmp <- gsub("Lymph node", "lymph_node",tmp)
tmp <- gsub("Soft Tissue", "soft_tissue",tmp)
tmp <- gsub("Bladder", "bladder",tmp)
tmp <- gsub("Liver", "liver",tmp)
tmp <- gsub("Prostate", "prostate",tmp)
tmp <- gsub("Psoas Mass", "pelvic",tmp)
tmp <- gsub("Retroperitonium", "peritoneal",tmp)
tmp <- gsub("Neck Mass", "neck",tmp)
tmp <- gsub("Dura", "brain",tmp)
tmp <- gsub("Abdominal wall mass", "peritoneal",tmp)
tmp <- gsub("Peritoneal nodule", "peritoneal",tmp)
tmp <- gsub("Perirectal", "rectal",tmp)
tmp <- gsub("Thoracic epidural", "thoracic",tmp)
tmp <- gsub("Chest Wall", "peritoneal",tmp)
tmp <- gsub("Penile", "penile",tmp)
tmp <- gsub("ProstateTURP", "prostate",tmp)
tmp <- gsub("Muscle", "muscle",tmp)
tmp <- gsub("glans penis", "penile",tmp)
tmp <- gsub("Pelvic mass", "pelvic",tmp)
curated$metatstatic_site <- tmp



#get the second file with more sample specfic information
uncurated <- read.delim("../DFCI_version/uncurated_copy/prad_su2c_2015/data_clinical_sample.txt",as.is=TRUE,row.names=1,sep="\t", skip=4)
colnames(uncurated) <- as.character(read.delim("../DFCI_version/uncurated_copy/prad_su2c_2015/data_clinical_sample.txt",as.is=TRUE,header=T,row.names=1,sep="\t",skip=3)[1,])
uncurated <- uncurated[row.names(curated),]

##tumor purity
curated$tumor_purity_pathology <- uncurated$TUMOR_CONTENT

write.table(curated, row.names = FALSE, file="../DFCI_version/curated/prad_su2c_2015_curated_pdata.txt",sep="\t")

