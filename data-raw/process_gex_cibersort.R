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

load("data-raw/gex_true.RData")
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




