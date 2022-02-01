##################################################################################################################
#Mutation - Raggedexp
##################################################################################################################

#TCGA

# Download the data and unzip
UCSCXenaTools::XenaGenerate(subset = XenaHostNames=="tcgaHub") %>% 
  UCSCXenaTools::XenaFilter(filterDatasets = "./data-raw/mc3/PRAD_mc3.txt") %>% 
  UCSCXenaTools::XenaFilter(filterDatasets = "PRAD") -> df_todo

UCSCXenaTools::XenaQuery(df_todo) %>%
  UCSCXenaTools::XenaDownload(destdir = "./data-raw") -> xe_download

R.utils::gunzip("./data-raw/mc3/PRAD_mc3.txt.gz")

# Process it and store it as raggedexperiment object
tcga_mut<-rio::import("./data-raw/mc3/PRAD_mc3.txt")

tcga_mut<-tcga_mut[,c(2:4,1,5:12)]
colnames(tcga_mut)[1:3]=c("seqnames","start","end")

tcga_mut$sample<-gsub("-",".",tcga_mut$sample)
names(tcga_mut)[names(tcga_mut) == 'effect'] <- "Variant_Classification"

a=subset(tcga_mut, sample %in% colnames(mae_tcga[["gex.tpm"]]))

GRL <- makeGRangesListFromDataFrame(a, split.field = "sample",
                                    names.field = "gene",keep.extra.columns = TRUE)

ragexp_tcga=RaggedExperiment::RaggedExperiment(GRL)


#mae_tcga[["mut"]]<-NULL
mae_tcga <- c(mae_tcga, mut = ragexp_tcga)
usethis::use_data(mae_tcga, overwrite = TRUE)


#Ren et al
ren_mut=generate_cbioportaldata_mut("prad_eururol_2017")
ragexp_ren=ren_mut[["mutations"]]

mae_ren <- c(mae_ren, mut = ragexp_ren)
usethis::use_data(mae_ren, overwrite = TRUE)

#ren_mut=generate_cbioportal_mut("ren",genes=curatedPCaData::curatedPCaData_genes$hgnc_symbol)


# ren_mut2<-rio::import("./data-raw/data_mutations_extended_ren.txt")
# 
# ren_mut2<-ren_mut2[,c(5:46,1:4)]
# 
# colnames(ren_mut2)[1:4]=c("seqnames","start","end","strand")
# 
# GRL <- makeGRangesListFromDataFrame(ren_mut2, split.field = "Tumor_Sample_Barcode",
#                                     names.field = "Hugo_Symbol",keep.extra.columns = TRUE)
# 
# 
# 
# ragexp_ren=RaggedExperiment::RaggedExperiment(GRL)


#mae_ren[["mut"]]<-NULL


#as(ragexp_ren, "GRangesList")

#Barbieri et al

barbieri_mut=generate_cbioportaldata_mut("prad_broad")
ragexp_barbieri=barbieri_mut[["mutations"]]
colnames(ragexp_barbieri)<-gsub("-",".",colnames(ragexp_barbieri))

save(ragexp_barbieri, file="data-raw/mut_barbieri.RData")
mae_barbieri <- c(mae_barbieri, mut = ragexp_barbieri)
usethis::use_data(mae_barbieri, overwrite = TRUE)

# barbieri_mut<-rio::import("./data-raw/data_mutations_extended_barbieri.txt")
# barbieri_mut2<-barbieri_mut[,c(5:113,1:4)]
# colnames(barbieri_mut2)[1:4]=c("seqnames","start","end","strand")
# barbieri_mut2$Tumor_Sample_Barcode<-gsub("-", ".", barbieri_mut2$Tumor_Sample_Barcode, fixed=TRUE)
# 
# GRL <- makeGRangesListFromDataFrame(barbieri_mut2, split.field = "Tumor_Sample_Barcode",
#                                     names.field = "Hugo_Symbol",keep.extra.columns = TRUE)
# 
# ragexp_barbieri=RaggedExperiment::RaggedExperiment(GRL)

#mae_barbieri[["mut"]]<-NULL


#Abida et al


abida_mut=curatedPCaData::generate_cbioportaldata_mut("prad_su2c_2019")
ragexp_abida=abida_mut[["mutations"]]
colnames(ragexp_abida)<-gsub("-",".",colnames(ragexp_abida))

mae_abida <- c(mae_abida, mut = ragexp_abida)
usethis::use_data(mae_abida, overwrite = TRUE)


# mycancerstudy = cgdsr::getCancerStudies(mycgds)
# mycaselistren = cgdsr::getCaseLists(mycgds,"prad_su2c_2019")
# gp = cgdsr::getGeneticProfiles(mycgds,"prad_su2c_2019")
# 
# 
# u<-cgdsr::getMutationData(mycgds,"prad_su2c_2019","prad_su2c_2019_mutations",c(d))
# abida_mut<-rio::import("./data-raw/data_mutations_extended_abida.txt")
# abida_mut2<-abida_mut[,c(5:46,1:4)]
# colnames(abida_mut2)[1:4]=c("seqnames","start","end","strand")
# abida_mut2$Tumor_Sample_Barcode<-gsub("-", ".", abida_mut2$Tumor_Sample_Barcode, fixed=TRUE)
# 
# GRL <- makeGRangesListFromDataFrame(abida_mut2, split.field = "Tumor_Sample_Barcode",
#                                     names.field = "Hugo_Symbol",keep.extra.columns = TRUE)
# 
# ragexp_abida=RaggedExperiment::RaggedExperiment(GRL)

#mae_abida[["mut"]]<-NULL


#Taylor et al

#mycgds <- cgdsr::CGDS("http://www.cbioportal.org/")
#prad_mskcc_mutations
#uncurated_cbio <- cgdsr::getMutationData.CGDS(mycgds, caseList = "prad_mskcc",geneticProfile="prad_mskcc_mutations",genes=taylor_mut$Hugo_Symbol)
taylor_mut<-curatedPCaData:::generate_cbioportaldata_mut("prad_mskcc")
ragexp_taylor<-taylor_mut[["mutations"]]
same_barcode=colnames(ragexp_taylor)[grepl("PCA", colnames(ragexp_taylor))]
ind=which(colnames(ragexp_taylor) %in% same_barcode=="TRUE")
ragexp_taylor2<-ragexp_taylor[, c(ind)]

mae_taylor <- c(mae_taylor, mut = ragexp_taylor2)
usethis::use_data(mae_taylor, overwrite = TRUE)


#colnames(ragexp_taylor)=colnames(grepl("PCA", colnames(ragexp_taylor)))
# t=ragexp_taylor
# taylor_mut<-rio::import("./data-raw/data_mutations_extended_taylor.txt")
# #taylor_mut2<-taylor_mut[,c(17,5:16,18:160,1:4)]
# #taylor_mut2<-taylor_mut[,c(5:160,1:4)]
# colnames(taylor_mut)[5:8]=c("seqnames","start","end","strand")
# 
# same_barcode=taylor_mut[grepl("PCA", taylor_mut$Tumor_Sample_Barcode),]
# same_barcode=same_barcode[,c(5:8,1,9,2:4,10:160)]
# 
# GRL <- makeGRangesListFromDataFrame(same_barcode, split.field = "Tumor_Sample_Barcode",
#                                     names.field = "Hugo_Symbol",keep.extra.columns = TRUE)
# ragexp_taylor=RaggedExperiment(GRL)
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

#mae_taylor[["mut"]]<-NULL

