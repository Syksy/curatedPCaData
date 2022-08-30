#############################################################################################################
# NEW CIBERSORT
#############################################################################################################
gex=mae_barwick[["gex.logq"]]
gex2 <- gex
gex2 <- as.data.frame(gex2)
gex2<- gex2[complete.cases(gex2), ]
#gex_chandran2 <- gex_chandran2[-9561,]
gex2[,147]<-rownames(gex2)
#View(gex2)
rownames(gex2) <- NULL
gex2 <- gex2[,c(147,1:146)]
#gex2 <- gex2[-9534,]
write.table(gex2,"/Users/varsha/Desktop/gex_cibersortx/changed_gex/gex_barwick.txt",sep="\t",row.names = FALSE,col.names = TRUE)


load("data-raw/gex.logr_kunderfranco.RData")
gex_kunderfranco2 <- gex.logr_kunderfranco
gex_kunderfranco2 <- as.data.frame(gex_kunderfranco2)
gex_kunderfranco2<- gex_kunderfranco2[complete.cases(gex_kunderfranco2), ]
#gex_kunderfranco_chandran2 <- gex_kunderfranco_chandran2[-9561,]
gex_kunderfranco2[,68]<-rownames(gex_kunderfranco2)
#View(gex_kunderfranco2)
rownames(gex_kunderfranco2) <- NULL
gex_kunderfranco2 <- gex_kunderfranco2[,c(68,1:67)]
#gex_kunderfranco2 <- gex_kunderfranco2[-9534,]
write.table(gex_kunderfranco2,"/Users/varsha/Desktop/gex_cibersortx/changed_gex/gex_kunderfranco.txt",sep="\t",row.names = FALSE,col.names = TRUE)


load("data-raw/gex.logr_true.RData")
gex_true2 <- gex
gex_true2 <- as.data.frame(gex_true2)
gex_true2<- gex_true2[complete.cases(gex_true2), ]
#gex_true2 <- gex_true2[-9561,]
gex_true2[,33]<-rownames(gex_true2)
#View(gex_true2)
rownames(gex_true2) <- NULL
gex_true2 <- gex_true2[,c(33,1:32)]
#gex_true2 <- gex_true2[-9534,]
write.table(gex_true2,"/Users/varsha/Desktop/gex_cibersortx/changed_gex/gex_true.txt",sep="\t",row.names = FALSE,col.names = TRUE)

gex_Barbieri<-curatedPCaData::mae_barbieri[["gex.relz"]]
gex_barbieri2 <- gex_Barbieri
gex_barbieri2 <- as.data.frame(gex_barbieri2)
gex_barbieri2<- gex_barbieri2[complete.cases(gex_barbieri2), ]
#gex_barbieri2 <- gex_barbieri2[-9561,]
gex_barbieri2[,32]<-rownames(gex_barbieri2)
#View(gex_barbieri2)
rownames(gex_barbieri2) <- NULL
gex_barbieri2 <- gex_barbieri2[,c(32,1:31)]
#gex_barbieri2 <- gex_barbieri2[-9534,]
write.table(gex_barbieri2,"/Users/varsha/Desktop/gex_cibersortx/changed_gex_latest/gex_barbieri.txt",sep="\t",row.names = FALSE,col.names = TRUE)

gex_eururol<-curatedPCaData::mae_ren[["gex.relz"]]
gex_ren2 <- gex_eururol
gex_ren2 <- as.data.frame(gex_ren2)
gex_ren2<- gex_ren2[complete.cases(gex_ren2), ]
#gex_ren2 <- gex_ren2[-9561,]
gex_ren2[,66]<-rownames(gex_ren2)
#View(gex_ren2)
rownames(gex_ren2) <- NULL
gex_ren2 <- gex_ren2[,c(66,1:65)]
#gex_ren2 <- gex_ren2[-9534,]
write.table(gex_ren2,"/Users/varsha/Desktop/gex_cibersortx/changed_gex_latest/gex_ren.txt",sep="\t",row.names = FALSE,col.names = TRUE)

load("data-raw/gex.rma_igc.RData")
gex_igc2 <- gex.rma_igc
gex_igc2 <- as.data.frame(gex_igc2)
gex_igc2<- gex_igc2[complete.cases(gex_igc2), ]
#gex_igc2 <- gex_igc2[-9561,]
gex_igc2[,84]<-rownames(gex_igc2)
#View(gex_igc2)
rownames(gex_igc2) <- NULL
gex_igc2 <- gex_igc2[,c(84,1:83)]
#gex_igc2 <- gex_igc2[-6860,]
write.table(gex_igc2,"/Users/varsha/Desktop/gex_cibersortx/changed_gex/gex_igc.txt",sep="\t",row.names = FALSE,col.names = TRUE)

load("data-raw/gex.rma_kim.RData")
gex_kim2 <- gex_kim
gex_kim2 <- as.data.frame(gex_kim2)
gex_kim2<- gex_kim2[complete.cases(gex_kim2), ]
#gex_kim2 <- gex_kim2[-9561,]
gex_kim2[,267]<-rownames(gex_kim2)
#View(gex_kim2)
rownames(gex_kim2) <- NULL
gex_kim2 <- gex_kim2[,c(267,1:266)]
#gex_kim2 <- gex_kim2[-6860,]
write.table(gex_kim2,"/Users/varsha/Desktop/gex_cibersortx/changed_gex/gex_kim.txt",sep="\t",row.names = FALSE,col.names = TRUE)

load("data-raw/gex.rma_taylor.RData")
gex_taylor2 <- gex
gex_taylor2 <- as.data.frame(gex_taylor2)
gex_taylor2<- gex_taylor2[complete.cases(gex_taylor2), ]
#gex_taylor2 <- gex_taylor2[-9561,]
gex_taylor2[,180]<-rownames(gex_taylor2)
#View(gex_taylor2)
rownames(gex_taylor2) <- NULL
gex_taylor2 <- gex_taylor2[,c(180,1:179)]
#gex_taylor2 <- gex_taylor2[-6860,]
write.table(gex_taylor2,"/Users/varsha/Desktop/gex_cibersortx/changed_gex/gex_taylor.txt",sep="\t",row.names = FALSE,col.names = TRUE)

gex.rma_wallace=mae_wallace[["gex.rma"]]
#gex.rma_wallace<-mae_wallace[["gex.rma"]]
gex_wallace2 <- gex.rma_wallace
gex_wallace2 <- as.data.frame(gex_wallace2)
gex_wallace2<- gex_wallace2[complete.cases(gex_wallace2), ]
#gex_wallace2 <- gex_wallace2[-9561,]
gex_wallace2[,84]<-rownames(gex_wallace2)
#View(gex_wallace2)
rownames(gex_wallace2) <- NULL
gex_wallace2 <- gex_wallace2[,c(84,1:83)]
#gex_wallace2 <- gex_wallace2[-6860,]
write.table(gex_wallace2,"/Users/varsha/Desktop/gex_cibersortx/changed_gex/gex_wallace.txt",sep="\t",row.names = FALSE,col.names = TRUE)

load("data-raw/gex.rma_sun.RData")
gex_sun2 <- gex
gex_sun2 <- as.data.frame(gex_sun2)
#gex_sun2 <- gex_sun2[-9561,]
gex_sun2[,80]<-rownames(gex_sun2)
#View(gex_sun2)
rownames(gex_sun2) <- NULL
gex_sun2 <- gex_sun2[,c(80,1:79)]
#gex_sun2 <- gex_sun2[-9534,]
write.table(gex_sun2,"/Users/varsha/Desktop/gex_cibersortx/changed_gex/gex_sun.txt",sep="\t",row.names = FALSE,col.names = TRUE)

load("data-raw/gex.rma_wang.RData")
gex_wang2 <- gex
gex_wang2 <- as.data.frame(gex_wang2)
gex_wang2<- gex_wang2[complete.cases(gex_wang2), ]
#gex_wang2 <- gex_wang2[-9561,]
gex_wang2[,149]<-rownames(gex_wang2)
#View(gex_wang2)
rownames(gex_wang2) <- NULL
gex_wang2 <- gex_wang2[,c(149,1:148)]
#gex_wang2 <- gex_wang2[-6860,]
write.table(gex_wang2,"/Users/varsha/Desktop/gex_cibersortx/changed_gex/gex_wang.txt",sep="\t",row.names = FALSE,col.names = TRUE)

load("data-raw/gex.logq_friedrich.RData")
gex_friedrich2 <- gex.logq_friedrich
gex_friedrich2 <- as.data.frame(gex_friedrich2)
gex_friedrich2<- gex_friedrich2[complete.cases(gex_friedrich2), ]
#gex_friedrich2 <- gex_friedrich2[-9561,]
gex_friedrich2[,256]<-rownames(gex_friedrich2)
#View(gex_friedrich2)
rownames(gex_friedrich2) <- NULL
gex_friedrich2 <- gex_friedrich2[,c(256,1:255)]
#gex_friedrich2 <- gex_friedrich2[-6860,]
write.table(gex_friedrich2,"/Users/varsha/Desktop/gex_cibersortx/changed_gex/gex_friedrich.txt",sep="\t",row.names = FALSE,col.names = TRUE)

gex_polyA_abida=curatedPCaData::mae_abida[["gex.relz"]]
gex_abida2 <- gex_polyA_abida
gex_abida2 <- as.data.frame(gex_abida2)
gex_abida2<- gex_abida2[complete.cases(gex_abida2), ]
#gex_abida2 <- gex_abida2[-9561,]
gex_abida2[,267]<-rownames(gex_abida2)
#View(gex_abida2)
rownames(gex_abida2) <- NULL
gex_abida2 <- gex_abida2[,c(267,1:266)]
#gex_abida2 <- gex_abida2[-6860,]
write.table(gex_abida2,"/Users/varsha/Desktop/gex_cibersortx/changed_gex_latest/gex_abida.txt",sep="\t",row.names = FALSE,col.names = TRUE)



load("data-raw/gex.rma_chandran.RData")
gex_chandran2 <- gex.rma_chandran
gex_chandran2 <- as.data.frame(gex_chandran2)
gex_chandran2<- gex_chandran2[complete.cases(gex_chandran2), ]
#gex_chandran2 <- gex_chandran2[-9561,]
gex_chandran2[,172]<-rownames(gex_chandran2)
#View(gex_chandran2)
rownames(gex_chandran2) <- NULL
gex_chandran2 <- gex_chandran2[,c(172,1:171)]
#gex_chandran2 <- gex_chandran2[-9534,]
write.table(gex_chandran2,"/Users/varsha/Desktop/gex_cibersortx/changed_gex/gex_chandran.txt",sep="\t",row.names = FALSE,col.names = TRUE)

load("data-raw/gex.rma_weiner.RData")
gex_weiner2 <- gex.rma_weiner
gex_weiner2 <- as.data.frame(gex_weiner2)
gex_weiner2<- gex_weiner2[complete.cases(gex_weiner2), ]
#gex_weiner2 <- gex_weiner2[-9561,]
gex_weiner2[,839]<-rownames(gex_weiner2)
#View(gex_weiner2)
rownames(gex_weiner2) <- NULL
gex_weiner2 <- gex_weiner2[,c(839,1:838)]
#gex_weiner2 <- gex_weiner2[-9534,]
write.table(gex_weiner2,"/Users/varsha/Desktop/gex_cibersortx/changed_gex/gex_weiner.txt",sep="\t",row.names = FALSE,col.names = TRUE)


load("data-raw/gex.rma_icgcca.RData")
gex_icgcca2 <- gex_icgcca
gex_icgcca2 <- as.data.frame(gex_icgcca2)
gex_icgcca2<- gex_icgcca2[complete.cases(gex_icgcca2), ]
#gex_icgcca2 <- gex_icgcca2[-9561,]
gex_icgcca2[,214]<-rownames(gex_icgcca2)
#View(gex_icgcca2)
rownames(gex_icgcca2) <- NULL
gex_icgcca2 <- gex_icgcca2[,c(214,1:213)]
#gex_icgcca2 <- gex_icgcca2[-6860,]
write.table(gex_icgcca2,"/Users/varsha/Desktop/gex_cibersortx/changed_gex/gex_icgcca.txt",sep="\t",row.names = FALSE,col.names = TRUE)

gex.rsem.log_tcga=curatedPCaData::mae_tcga[["gex.rsem.log"]]
gex_tcga2 <- gex.rsem.log_tcga
gex_tcga2 <- as.data.frame(gex_tcga2)
gex_tcga2<- gex_tcga2[complete.cases(gex_tcga2), ]
#gex_tcga2 <- gex_tcga2[-9561,]
gex_tcga2[,334]<-rownames(gex_tcga2)
#View(gex_tcga2)
rownames(gex_tcga2) <- NULL
gex_tcga2 <- gex_tcga2[,c(334,1:333)]
#gex_tcga2 <- gex_tcga2[-9534,]
write.table(gex_tcga2,"/Users/varsha/Desktop/gex_cibersortx/changed_gex_latest/gex_tcga.txt",sep="\t",row.names = FALSE,col.names = TRUE)




