##################################################################################################################
#Mutation - Raggedexp
##################################################################################################################

#TCGA

# Download the data and unzip
XenaGenerate(subset = XenaHostNames=="tcgaHub") %>% 
  XenaFilter(filterDatasets = "mc3/PRAD_mc3.txt") %>% 
  XenaFilter(filterDatasets = "PRAD") -> df_todo

XenaQuery(df_todo) %>%
  XenaDownload(destdir = "./data-raw") -> xe_download

R.utils::gunzip("./data-raw/mc3/PRAD_mc3.txt.gz")

# Process it and store it as raggedexperiment object
tcga_mut<-rio::import("./data-raw/mc3/PRAD_mc3.txt")

tcga_mut<-tcga_mut[,c(2:4,1,5:12)]
colnames(tcga_mut)[1:3]=c("seqnames","start","end")

tcga_mut$sample<-gsub("-",".",tcga_mut$sample)
names(tcga_mut)[names(tcga_mut) == 'effect'] <- "Variant_Classification"

GRL <- makeGRangesListFromDataFrame(tcga_mut, split.field = "sample",
                                    names.field = "gene",keep.extra.columns = TRUE)

ragexp_tcga=RaggedExperiment::RaggedExperiment(GRL)


mae_tcga[["mut_ragex"]]<-NULL
mae_tcga <- c(mae_tcga, mut_ragex = ragexp_tcga)
usethis::use_data(mae_tcga, overwrite = TRUE)


#Ren et al

ren_mut<-rio::import("/Users/varsha/Downloads/prad_eururol_2017/data_mutations_extended.txt")
#ren_mut=cgdsr::getMutationData(mycgds, caseList="prad_eururol_2017_sequenced")
ren_mut2<-ren_mut[,c(5:46,1:4)]
colnames(ren_mut2)[1:4]=c("seqnames","start","end","strand")

GRL <- makeGRangesListFromDataFrame(ren_mut2, split.field = "Tumor_Sample_Barcode",
                                    names.field = "Hugo_Symbol",keep.extra.columns = TRUE)

ragexp_ren=RaggedExperiment::RaggedExperiment(GRL)


mae_ren[["mut_ragex"]]<-NULL
mae_ren <- c(mae_ren, mut_ragex = ragexp_ren)
usethis::use_data(mae_ren, overwrite = TRUE)

#as(ragexp_ren, "GRangesList")

#Barbieri et al

barbieri_mut<-rio::import("/Users/varsha/Downloads/prad_broad/data_mutations_extended.txt")
barbieri_mut2<-barbieri_mut[,c(5:113,1:4)]
colnames(barbieri_mut2)[1:4]=c("seqnames","start","end","strand")
barbieri_mut2$Tumor_Sample_Barcode<-gsub("-", ".", barbieri_mut2$Tumor_Sample_Barcode, fixed=TRUE)

GRL <- makeGRangesListFromDataFrame(barbieri_mut2, split.field = "Tumor_Sample_Barcode",
                                    names.field = "Hugo_Symbol",keep.extra.columns = TRUE)

ragexp_barbieri=RaggedExperiment::RaggedExperiment(GRL)

mae_barbieri[["mut_ragex"]]<-NULL
mae_barbieri <- c(mae_barbieri, mut_ragex = ragexp_barbieri)
usethis::use_data(mae_barbieri, overwrite = TRUE)

#Abida et al

abida_mut<-rio::import("/Users/varsha/Downloads/prad_su2c_2019/data_mutations_extended.txt")
# mycancerstudy = cgdsr::getCancerStudies(mycgds)
# mycaselistren = cgdsr::getCaseLists(mycgds,"prad_su2c_2019")
# gp = cgdsr::getGeneticProfiles(mycgds,"prad_su2c_2019")
# 
# 
# u<-cgdsr::getMutationData(mycgds,"prad_su2c_2019","prad_su2c_2019_mutations",c(d))

abida_mut2<-abida_mut[,c(5:46,1:4)]
colnames(abida_mut2)[1:4]=c("seqnames","start","end","strand")
abida_mut2$Tumor_Sample_Barcode<-gsub("-", ".", abida_mut2$Tumor_Sample_Barcode, fixed=TRUE)

GRL <- makeGRangesListFromDataFrame(abida_mut2, split.field = "Tumor_Sample_Barcode",
                                    names.field = "Hugo_Symbol",keep.extra.columns = TRUE)

ragexp_abida=RaggedExperiment::RaggedExperiment(GRL)

mae_abida[["mut_ragex"]]<-NULL
mae_abida <- c(mae_abida, mut_ragex = ragexp_abida)
usethis::use_data(mae_abida, overwrite = TRUE)

#Taylor et al

mycgds <- cgdsr::CGDS("http://www.cbioportal.org/")
#prad_mskcc_mutations
uncurated_cbio <- cgdsr::getMutationData(mycgds, caseList = "prad_mskcc",geneticProfile="prad_mskcc_mutations",genes=c("FGF1","JAK2"))

taylor_mut<-rio::import("/Users/varsha/Downloads/prad_mskcc/data_mutations_extended.txt")
#taylor_mut2<-taylor_mut[,c(17,5:16,18:160,1:4)]
#taylor_mut2<-taylor_mut[,c(5:160,1:4)]
colnames(taylor_mut)[5:8]=c("seqnames","start","end","strand")

same_barcode=taylor_mut[grepl("PCA", taylor_mut$Tumor_Sample_Barcode),]
same_barcode=same_barcode[,c(5:8,1,9,2:4,10:160)]

GRL <- makeGRangesListFromDataFrame(same_barcode, split.field = "Tumor_Sample_Barcode",
                                    names.field = "Hugo_Symbol",keep.extra.columns = TRUE)
ragexp_taylor=RaggedExperiment(GRL)
#X<-split(same_barcode, same_barcode$Tumor_Sample_Barcode)

#Y<-split(same_barcode, same_barcode$Hugo_Symbol)
# count=1:92
# for (i in count){
# grl <- GRangesList(same_barcode$Hugo_Symbol[i] = X[[i]])
# }
#X$gene=as.list(same_barcode$Hugo_Symbol)
# count=1:43
# 
# for (i in count){
#   rownames(X[[i]])=make.names(X[[i]]$Hugo_Symbol, unique=TRUE)
# }
# 
# 
# b=GRangesList(X)
# ragexp_taylor2=RaggedExperiment::RaggedExperiment(b)

mae_taylor[["mut_ragex"]]<-NULL
mae_taylor <- c(mae_taylor, mut_ragex = ragexp_taylor)
usethis::use_data(mae_taylor, overwrite = TRUE)
