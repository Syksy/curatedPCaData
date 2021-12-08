#############################################################################################################
# NEW CIBERSORT
#############################################################################################################
load("data-raw/gex.logq_barwick.RData")
gex2 <- gex
gex2 <- as.data.frame(gex2)
gex2<- gex2[complete.cases(gex2), ]
#gex_chandran2 <- gex_chandran2[-9561,]
gex2[,181]<-rownames(gex2)
#View(gex2)
rownames(gex2) <- NULL
gex2 <- gex2[,c(181,1:180)]
#gex2 <- gex2[-9534,]
write.table(gex2,"/Users/varsha/Desktop/gex_cibersortx/changed_gex/gex_barwick.txt",sep="\t",row.names = FALSE,col.names = TRUE)


load("data-raw/gex_kunderfranco.RData")
gex_kunderfranco2 <- gex_kunderfranco
gex_kunderfranco2 <- as.data.frame(gex_kunderfranco2)
gex_kunderfranco2<- gex_kunderfranco2[complete.cases(gex_kunderfranco2), ]
#gex_kunderfranco_chandran2 <- gex_kunderfranco_chandran2[-9561,]
gex_kunderfranco2[,68]<-rownames(gex_kunderfranco2)
#View(gex_kunderfranco2)
rownames(gex_kunderfranco2) <- NULL
gex_kunderfranco2 <- gex_kunderfranco2[,c(68,1:67)]
#gex_kunderfranco2 <- gex_kunderfranco2[-9534,]
write.table(gex_kunderfranco2,"/Users/varsha/Desktop/gex_cibersortx/changed_gex/gex_kunderfranco.txt",sep="\t",row.names = FALSE,col.names = TRUE)


load("data-raw/gex_true.RData")
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

load("data-raw/gex.relz_barbieri.RData")
gex_barbieri2 <- gex_Barbieri
gex_barbieri2 <- as.data.frame(gex_barbieri2)
gex_barbieri2<- gex_barbieri2[complete.cases(gex_barbieri2), ]
#gex_barbieri2 <- gex_barbieri2[-9561,]
gex_barbieri2[,21]<-rownames(gex_barbieri2)
#View(gex_barbieri2)
rownames(gex_barbieri2) <- NULL
gex_barbieri2 <- gex_barbieri2[,c(21,1:20)]
#gex_barbieri2 <- gex_barbieri2[-9534,]
write.table(gex_barbieri2,"/Users/varsha/Desktop/gex_cibersortx/changed_gex/gex_barbieri.txt",sep="\t",row.names = FALSE,col.names = TRUE)

load("data-raw/gex.relz_ren.RData")
gex_ren2 <- gex_eururol
gex_ren2 <- as.data.frame(gex_ren2)
gex_ren2<- gex_ren2[complete.cases(gex_ren2), ]
#gex_ren2 <- gex_ren2[-9561,]
gex_ren2[,66]<-rownames(gex_ren2)
#View(gex_ren2)
rownames(gex_ren2) <- NULL
gex_ren2 <- gex_ren2[,c(66,1:65)]
#gex_ren2 <- gex_ren2[-9534,]
write.table(gex_ren2,"/Users/varsha/Desktop/gex_cibersortx/changed_gex/gex_ren.txt",sep="\t",row.names = FALSE,col.names = TRUE)

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

load("data-raw/gex.rma_wallace.RData")
gex.rma_wallace<-mae_wallace[["gex.rma"]]
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

load("data-raw/gex.relz_abida.RData")
gex_abida2 <- gex_polyA_abida
gex_abida2 <- as.data.frame(gex_abida2)
gex_abida2<- gex_abida2[complete.cases(gex_abida2), ]
#gex_abida2 <- gex_abida2[-9561,]
gex_abida2[,267]<-rownames(gex_abida2)
#View(gex_abida2)
rownames(gex_abida2) <- NULL
gex_abida2 <- gex_abida2[,c(267,1:266)]
#gex_abida2 <- gex_abida2[-6860,]
write.table(gex_abida2,"/Users/varsha/Desktop/gex_cibersortx/changed_gex/gex_abida.txt",sep="\t",row.names = FALSE,col.names = TRUE)



######################################################################################################################

load("data-raw/gex_sun.RData")
gex_sun2 <- gex_sun
gex_sun2 <- as.data.frame(gex_sun2)
#gex_sun2 <- gex_sun2[-9561,]
gex_sun2[,80]<-rownames(gex_sun2)
#View(gex_sun2)
rownames(gex_sun2) <- NULL
gex_sun2 <- gex_sun2[,c(80,1:79)]
gex_sun2 <- gex_sun2[-9534,]
write.table(gex_sun2,"gex_sun.txt",sep="\t",row.names = FALSE,col.names = TRUE)

load("data-raw/gex_chandran.RData")
gex_chandran2 <- gex_chandran
gex_chandran2 <- as.data.frame(gex_chandran2)
gex_chandran2<- gex_chandran2[complete.cases(gex_chandran2), ]
#gex_chandran2 <- gex_chandran2[-9561,]
gex_chandran2[,504]<-rownames(gex_chandran2)
#View(gex_chandran2)
rownames(gex_chandran2) <- NULL
gex_chandran2 <- gex_chandran2[,c(504,1:503)]
#gex_chandran2 <- gex_chandran2[-9534,]
write.table(gex_chandran2,"gex_chandran.txt",sep="\t",row.names = FALSE,col.names = TRUE)

load("data-raw/gex_weiner.RData")
gex_weiner2 <- gex_weiner
gex_weiner2 <- as.data.frame(gex_weiner2)
#gex_weiner2 <- gex_weiner2[-9561,]
gex_weiner2[,839]<-rownames(gex_weiner2)
#View(gex_weiner2)
rownames(gex_weiner2) <- NULL
gex_weiner2 <- gex_weiner2[,c(839,1:838)]
#gex_weiner2 <- gex_weiner2[-9534,]
write.table(gex_weiner2,"gex_weiner.txt",sep="\t",row.names = FALSE,col.names = TRUE)

load("data-raw/gex.logr_true.RData")
gex_true2 <- gex_true
#gex_true2 <- as.data.frame(gex_true2)
gex_true2<- gex_true2[complete.cases(gex_true2), ]
#gex_true2 <- gex_true2[-9561,]
gex_true2[,33]<-rownames(gex_true2)
#View(gex_true2)
rownames(gex_true2) <- NULL
gex_true2 <- gex_true2[,c(33,1:32)]
#gex_true2 <- gex_true2[-9534,]
write.table(gex_true2,"gex_true.txt",sep="\t",row.names = FALSE,col.names = TRUE)

load("data-raw/gex_igc.RData")
gex_igc2 <- gex_igc
gex_igc2 <- as.data.frame(gex_igc2)
#gex_igc2 <- gex_igc2[-9561,]
gex_igc2[,84]<-rownames(gex_igc2)
#View(gex_igc2)
rownames(gex_igc2) <- NULL
gex_igc2 <- gex_igc2[,c(84,1:83)]
gex_igc2 <- gex_igc2[-6860,]
write.table(gex_igc2,"gex_igc.txt",sep="\t",row.names = FALSE,col.names = TRUE)

load("data-raw/gex_capture_abida.RData")
gex_capture_abida2 <- gex_capture_abida
gex_capture_abida2 <- as.data.frame(gex_capture_abida2)
#gex_capture_abida2 <- gex_capture_abida2[-9561,]
gex_capture_abida2[,209]<-rownames(gex_capture_abida2)
#View(gex_capture_abida2)
rownames(gex_capture_abida2) <- NULL
gex_capture_abida2 <- gex_capture_abida2[,c(209,1:208)]
#gex_capture_abida2 <- gex_capture_abida2[-6860,]
write.table(gex_capture_abida2,"gex_capture_abida.txt",sep="\t",row.names = FALSE,col.names = TRUE)

load("data-raw/gex.relz_abida.RData")
gex_polyA_abida2 <- gex_polyA_abida
gex_polyA_abida2 <- as.data.frame(gex_polyA_abida2)
#gex_polyA_abida2 <- gex_polyA_abida2[-9561,]
gex_polyA_abida2[,267]<-rownames(gex_polyA_abida2)
#View(gex_polyA_abida2)
rownames(gex_polyA_abida2) <- NULL
gex_polyA_abida2 <- gex_polyA_abida2[,c(267,1:266)]
#gex_polyA_abida2 <- gex_polyA_abida2[-6860,]
write.table(gex_polyA_abida2,"/Users/varsha/Desktop/gex_cibersortx/changed_gex/gex_polyA_abida.txt",sep="\t",row.names = FALSE,col.names = TRUE)

load("data-raw/gex_icgcca.RData")
gex_icgcca2 <- gex_icgcca
gex_icgcca2 <- as.data.frame(gex_icgcca2)
#gex_icgcca2 <- gex_icgcca2[-9561,]
gex_icgcca2[,214]<-rownames(gex_icgcca2)
#View(gex_icgcca2)
rownames(gex_icgcca2) <- NULL
gex_icgcca2 <- gex_icgcca2[,c(214,1:213)]
#gex_icgcca2 <- gex_icgcca2[-6860,]
write.table(gex_icgcca2,"gex_icgcca.txt",sep="\t",row.names = FALSE,col.names = TRUE)

load("data-raw/gex.logq_friedrich.RData")
gex_friedrich2 <- gex.logq_friedrich
gex_friedrich2 <- as.data.frame(gex_friedrich2)
#gex_friedrich2 <- gex_friedrich2[-9561,]
gex_friedrich2[,256]<-rownames(gex_friedrich2)
#View(gex_friedrich2)
rownames(gex_friedrich2) <- NULL
gex_friedrich2 <- gex_friedrich2[,c(256,1:255)]
#gex_friedrich2 <- gex_friedrich2[-6860,]
write.table(gex_friedrich2,"/Users/varsha/Desktop/gex_cibersortx/changed_gex/gex_friedrich.txt",sep="\t",row.names = FALSE,col.names = TRUE)




