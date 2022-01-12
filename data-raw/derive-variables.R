#####
##
## -- DERIVING NEW VARIABLES --
## Set of scripts that derive new variables from the original raw data
##
#####

## IMMUNE DECONVOLUTION VARIABLES

library(immunedeconv) # NECESSARY!


#####################################################################
#####################################################################
##                                                                 ##
##  The script is grouped by deconvolution method not by dataset   ##
##                                                                 ##   
#####################################################################
#####################################################################


## IMMUNE DECONVOLUTION VARIABLES

#####################################################
#####################################################
##                                                 ##
##                     CIBERSORTX                  ##
##                                                 ##
#####################################################
#####################################################

# Kunderfranco et al.

cibersort_kunderfranco<-rio::import("data-raw/CIBERSORTx_kunderfranco_Results.csv")
cibersort_kunderfranco<-cibersort_kunderfranco[ , -which(names(cibersort_kunderfranco) %in% c("P-value","Correlation","RMSE"))]
cibersort_kunderfranco <- t(cibersort_kunderfranco)
colnames(cibersort_kunderfranco)<- cibersort_kunderfranco[1,]
cibersort_kunderfranco<-cibersort_kunderfranco[-1,]

cibersort_kunderfranco<-as.matrix(cibersort_kunderfranco)
#save(cibersort_kunderfranco, file="data-raw/cibersort_kunderfranco.RData")
mae_kunderfranco <- c(mae_kunderfranco, cibersort = cibersort_kunderfranco)
#mae_kunderfranco <- create_mae(study_name = "kunderfranco")
usethis::use_data(mae_kunderfranco, overwrite = TRUE)

# Sun et al.

cibersort_sun<-rio::import("data-raw/CIBERSORTx_sun_Results.csv")
cibersort_sun<-cibersort_sun[ , -which(names(cibersort_sun) %in% c("P-value","Correlation","RMSE"))]
cibersort_sun <- t(cibersort_sun)
colnames(cibersort_sun)<- cibersort_sun[1,]
cibersort_sun<-cibersort_sun[-1,]
cibersort_sun<-as.matrix(cibersort_sun)
#save(cibersort_sun, file="data-raw/cibersort_sun.RData")
mae_sun <- c(mae_sun, cibersort = cibersort_sun)
#mae_sun <- create_mae(study_name = "sun")
usethis::use_data(mae_sun, overwrite = TRUE)

# Icgcca et al.

cibersort_icgcca<-rio::import("data-raw/CIBERSORTx_icgcca_Results.csv")
cibersort_icgcca<-cibersort_icgcca[ , -which(names(cibersort_icgcca) %in% c("P-value","Correlation","RMSE"))]
cibersort_icgcca <- t(cibersort_icgcca)
colnames(cibersort_icgcca)<- cibersort_icgcca[1,]
cibersort_icgcca<-cibersort_icgcca[-1,]
cibersort_icgcca<-as.matrix(cibersort_icgcca)
#save(cibersort_icgcca, file="data-raw/cibersort_icgcca.RData")
#mae_icgcca[["cibersort"]]<-NULL
mae_icgcca <- c(mae_icgcca, cibersort = cibersort_icgcca)
#mae_icgcca <- create_mae(study_name = "icgcca")
usethis::use_data(mae_icgcca, overwrite = TRUE)

# Abida et al.


cibersort_polyA_abida<-rio::import("data-raw/CIBERSORTx_abida_Results.csv")
cibersort_polyA_abida<-cibersort_polyA_abida[ , -which(names(cibersort_polyA_abida) %in% c("P-value","Correlation","RMSE"))]
cibersort_polyA_abida <- t(cibersort_polyA_abida)
colnames(cibersort_polyA_abida)<- cibersort_polyA_abida[1,]
cibersort_polyA_abida<-cibersort_polyA_abida[-1,]
cibersort_polyA_abida<-as.matrix(cibersort_polyA_abida)
#save(cibersort_polyA_abida, file="data-raw/cibersort_polyA_abida.RData")
mae_abida <- c(mae_abida, cibersort = cibersort_polyA_abida)
#mae_capture_abida <- create_mae(study_name = "capture_abida")
usethis::use_data(mae_abida, overwrite = TRUE)

# Wang et al.

cibersort_wang<-rio::import("data-raw/CIBERSORTx_wang_Results.csv")
cibersort_wang<-cibersort_wang[ , -which(names(cibersort_wang) %in% c("P-value","Correlation","RMSE"))]
cibersort_wang <- t(cibersort_wang)
colnames(cibersort_wang)<- cibersort_wang[1,]
cibersort_wang<-cibersort_wang[-1,]
cibersort_wang<-as.matrix(cibersort_wang)
#save(cibersort_wang, file="data-raw/cibersort_wang.RData")
mae_wang <- c(mae_wang, cibersort = cibersort_wang)
#mae_wang <- create_mae(study_name = "wang")
usethis::use_data(mae_wang, overwrite = TRUE)

# Kim et al.

cibersort_kim<-rio::import("data-raw/CIBERSORTx_kim_Results.csv")
cibersort_kim<-cibersort_kim[ , -which(names(cibersort_kim) %in% c("P-value","Correlation","RMSE"))]
cibersort_kim <- t(cibersort_kim)
colnames(cibersort_kim)<- cibersort_kim[1,]
cibersort_kim<-cibersort_kim[-1,]
cibersort_kim<-as.matrix(cibersort_kim)
#save(cibersort_kim, file="data-raw/cibersort_kim.RData")
mae_kim <- c(mae_kim, cibersort = cibersort_kim)
#mae_kim <- create_mae(study_name = "kim")
usethis::use_data(mae_kim, overwrite = TRUE)

# Barbieri et al.

cibersort_barbieri<-rio::import("data-raw/CIBERSORTx_barbieri_Results.csv")
cibersort_barbieri<-cibersort_barbieri[ , -which(names(cibersort_barbieri) %in% c("P-value","Correlation","RMSE"))]
cibersort_barbieri <- t(cibersort_barbieri)
colnames(cibersort_barbieri)<- cibersort_barbieri[1,]
cibersort_barbieri<-cibersort_barbieri[-1,]
cibersort_barbieri<-as.matrix(cibersort_barbieri)
#save(cibersort_barbieri, file="data-raw/cibersort_barbieri.RData")
mae_barbieri <- c(mae_barbieri, cibersort = cibersort_barbieri)
#mae_barbieri <- create_mae(study_name = "barbieri")
usethis::use_data(mae_barbieri, overwrite = TRUE)

# Ren et al.

cibersort_ren<-rio::import("data-raw/CIBERSORTx_ren_Results.csv")
cibersort_ren<-cibersort_ren[ , -which(names(cibersort_ren) %in% c("P-value","Correlation","RMSE"))]
cibersort_ren <- t(cibersort_ren)
colnames(cibersort_ren)<- cibersort_ren[1,]
cibersort_ren<-cibersort_ren[-1,]
cibersort_ren<-as.matrix(cibersort_ren)
#save(cibersort_ren, file="data-raw/cibersort_ren.RData")
mae_ren <- c(mae_ren, cibersort = cibersort_ren)
#mae_ren <- create_mae(study_name = "ren")
usethis::use_data(mae_ren, overwrite = TRUE)

# Wallace et al.

cibersort_wallace<-rio::import("data-raw/CIBERSORTx_wallace_Results.csv")
cibersort_wallace<-cibersort_wallace[ , -which(names(cibersort_wallace) %in% c("P-value","Correlation","RMSE"))]
cibersort_wallace <- t(cibersort_wallace)
colnames(cibersort_wallace)<- cibersort_wallace[1,]
cibersort_wallace<-cibersort_wallace[-1,]
cibersort_wallace<-as.matrix(cibersort_wallace)
#save(cibersort_wallace, file="data-raw/cibersort_wallace.RData")
mae_wallace <- c(mae_wallace, cibersort = cibersort_wallace)
#mae_wallace <- create_mae(study_name = "wallace")
usethis::use_data(mae_wallace, overwrite = TRUE)

# igc et al.

cibersort_igc<-rio::import("data-raw/CIBERSORTx_igc_Results.csv")
cibersort_igc<-cibersort_igc[ , -which(names(cibersort_igc) %in% c("P-value","Correlation","RMSE"))]
cibersort_igc <- t(cibersort_igc)
colnames(cibersort_igc)<- cibersort_igc[1,]
cibersort_igc<-cibersort_igc[-1,]
cibersort_igc<-as.matrix(cibersort_igc)
#save(cibersort_igc, file="data-raw/cibersort_igc.RData")
mae_igc <- c(mae_igc, cibersort = cibersort_igc)
#mae_igc <- create_mae(study_name = "igc")
usethis::use_data(mae_igc, overwrite = TRUE)

# TCGA et al.

cibersort_tcga<-rio::import("data-raw/CIBERSORTx_tcga_Results.csv")
cibersort_tcga<-cibersort_tcga[ , -which(names(cibersort_tcga) %in% c("P-value","Correlation","RMSE"))]
cibersort_tcga <- t(cibersort_tcga)
colnames(cibersort_tcga)<- cibersort_tcga[1,]
cibersort_tcga<-cibersort_tcga[-1,]
cibersort_tcga<-as.matrix(cibersort_tcga)
#save(cibersort_tcga, file="data-raw/cibersort_tcga.RData")
mae_tcga <- c(mae_tcga, cibersort = cibersort_tcga)
#mae_tcga <- create_mae(study_name = "tcga")
usethis::use_data(mae_tcga, overwrite = TRUE)

# Taylor et al.

cibersort_taylor<-rio::import("data-raw/CIBERSORTx_taylor_Results.csv")
cibersort_taylor<-cibersort_taylor[ , -which(names(cibersort_taylor) %in% c("P-value","Correlation","RMSE"))]
cibersort_taylor <- t(cibersort_taylor)
colnames(cibersort_taylor)<- cibersort_taylor[1,]
cibersort_taylor<-cibersort_taylor[-1,]
cibersort_taylor<-as.matrix(cibersort_taylor)
#save(cibersort_taylor, file="data-raw/cibersort_taylor.RData")
mae_taylor <- c(mae_taylor, cibersort = cibersort_taylor)
#mae_taylor <- create_mae(study_name = "taylor")
usethis::use_data(mae_taylor, overwrite = TRUE)

# True et al.

cibersort_true<-rio::import("data-raw/CIBERSORTx_true_Results.csv")
cibersort_true<-cibersort_true[ , -which(names(cibersort_true) %in% c("P-value","Correlation","RMSE"))]
cibersort_true <- t(cibersort_true)
colnames(cibersort_true)<- cibersort_true[1,]
cibersort_true<-cibersort_true[-1,]
cibersort_true<-as.matrix(cibersort_true)
#save(cibersort_true, file="data-raw/cibersort_true.RData")
mae_true <- c(mae_true, cibersort = cibersort_true)
#mae_true <- create_mae(study_name = "true")
usethis::use_data(mae_true, overwrite = TRUE)

# Barwick et al.

cibersort_barwick<-rio::import("data-raw/CIBERSORTx_barwick_Results.csv")
cibersort_barwick<-cibersort_barwick[ , -which(names(cibersort_barwick) %in% c("P-value","Correlation","RMSE"))]
cibersort_barwick <- t(cibersort_barwick)
colnames(cibersort_barwick)<- cibersort_barwick[1,]
cibersort_barwick<-cibersort_barwick[-1,]
cibersort_barwick<-as.matrix(cibersort_barwick)
#save(cibersort_barwick, file="data-raw/cibersort_barwick.RData")
mae_barwick <- c(mae_barwick, cibersort = cibersort_barwick)
#mae_barwick <- create_mae(study_name = "barwick")
usethis::use_data(mae_barwick, overwrite = TRUE)


# Chandran et al.

cibersort_chandran<-rio::import("data-raw/CIBERSORTx_chandran_Results.csv")
cibersort_chandran<-cibersort_chandran[ , -which(names(cibersort_chandran) %in% c("P-value","Correlation","RMSE"))]
cibersort_chandran <- t(cibersort_chandran)
colnames(cibersort_chandran)<- cibersort_chandran[1,]
cibersort_chandran<-cibersort_chandran[-1,]
cibersort_chandran<-as.matrix(cibersort_chandran)
#save(cibersort_chandran, file="data-raw/cibersort_chandran.RData")
mae_chandran <- c(mae_chandran, cibersort = cibersort_chandran)
#mae_chandran <- create_mae(study_name = "chandran")
usethis::use_data(mae_chandran, overwrite = TRUE)

# Friedrich et al.

cibersort_friedrich<-rio::import("data-raw/CIBERSORTx_friedrich_Results.csv")
cibersort_friedrich<-cibersort_friedrich[ , -which(names(cibersort_friedrich) %in% c("P-value","Correlation","RMSE"))]
cibersort_friedrich <- t(cibersort_friedrich)
colnames(cibersort_friedrich)<- cibersort_friedrich[1,]
cibersort_friedrich<-cibersort_friedrich[-1,]
cibersort_friedrich<-as.matrix(cibersort_friedrich)
#save(cibersort_friedrich, file="data-raw/cibersort_friedrich.RData")
mae_friedrich <- c(mae_friedrich, cibersort = cibersort_friedrich)
#mae_friedrich <- create_mae(study_name = "friedrich")
usethis::use_data(mae_friedrich, overwrite = TRUE)



library(immunedeconv) # NECESSARY!

#####################################################
#####################################################
##                                                 ##
##                     XCELL                       ##
##                                                 ##
#####################################################
#####################################################


# Abida et al.


tmp <- as.data.frame(immunedeconv::deconvolute(mae_abida[["gex.relz"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_abida <- c(mae_abida, xcell = tmp)

# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_abida, overwrite = TRUE)

# Barbieri et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_barbieri[["gex.relz"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_barbieri <- c(mae_barbieri, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_barbieri, overwrite = TRUE)

# Barwick et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_barwick[["gex.logq"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_barwick <- c(mae_barwick, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_barwick, overwrite = TRUE)

# Chandran et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_chandran[["gex.rma"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_chandran <- c(mae_chandran, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_chandran, overwrite = TRUE)

# Friedrich et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_friedrich[["gex.logq"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_friedrich <- c(mae_friedrich, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_friedrich, overwrite = TRUE)

# ICGCCA

tmp <- as.data.frame(immunedeconv::deconvolute(mae_icgcca[["gex.rma"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_icgcca <- c(mae_icgcca, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_icgcca, overwrite = TRUE)

# Kunderfranco et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_kunderfranco[["gex.logr"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_kunderfranco <- c(mae_kunderfranco, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_kunderfranco, overwrite = TRUE)

# Ren et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_ren[["gex.relz"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_ren <- c(mae_ren, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_ren, overwrite = TRUE)

# True et al.

# tmp <- as.data.frame(immunedeconv::deconvolute(mae_true[["gex"]], method="xcell"))
# rownames(tmp) <- tmp$cell_type
# # Omit cell type column and store only data of cell type populations
# tmp <- as.matrix(tmp[,-1])
# mae_true <- c(mae_true, xcell = tmp)
# # Save the derived new 'assay' types to the mae-object
# usethis::use_data(mae_true, overwrite = TRUE)

# Sun et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_sun[["gex.rma"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_sun <- c(curatedPCaData::mae_sun, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_sun, overwrite = TRUE)

# Taylor et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_taylor[["gex.rma"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_taylor <- c(curatedPCaData::mae_taylor, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_taylor, overwrite = TRUE)

# TCGA

tmp <- as.data.frame(immunedeconv::deconvolute(mae_tcga[["gex.tpm"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_tcga <- c(curatedPCaData::mae_tcga, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_tcga, overwrite = TRUE)



# Wang et al.

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_wang[["gex.rma"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_wang <- c(curatedPCaData::mae_wang, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_wang, overwrite = TRUE)

# Kim et al.

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_kim[["gex.rma"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_kim <- c(curatedPCaData::mae_kim, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_kim, overwrite = TRUE)

# Wallace et al.

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_wallace[["gex.rma"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_wallace <- c(mae_wallace, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_wallace, overwrite = TRUE)

# Igc 

tmp <- as.data.frame(immunedeconv::deconvolute(mae_igc[["gex.rma"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_igc <- c(mae_igc, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_igc, overwrite = TRUE)

# Weiner et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_weiner[["gex.rma"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_weiner <- c(mae_weiner, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_weiner, overwrite = TRUE)

# True et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_true[["gex.logr"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_true <- c(mae_true, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_true, overwrite = TRUE)

#####################################################
#####################################################
##                                                 ##
##                     EPIC                        ##
##                                                 ##
#####################################################
#####################################################

# Abida et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_abida[["gex.relz"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_abida <- c(mae_abida, epic = tmp)

# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_abida, overwrite = TRUE)

# Barbieri et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_barbieri[["gex.relz"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_barbieri <- c(mae_barbieri, epic = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_barbieri, overwrite = TRUE,internal = TRUE)

# Barwick et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_barwick[["gex.logq"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_barwick <- c(mae_barwick, epic = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_barwick, overwrite = TRUE)

# Chandran et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_chandran[["gex.rma"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_chandran <- c(mae_chandran, epic = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_chandran, overwrite = TRUE)

# Friedrich et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_friedrich[["gex.logq"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_friedrich <- c(mae_friedrich, epic = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_friedrich, overwrite = TRUE)

# ICGCCA

tmp <- as.data.frame(immunedeconv::deconvolute(mae_icgcca[["gex.rma"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_icgcca <- c(mae_icgcca, epic = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_icgcca, overwrite = TRUE)

# Kunderfranco et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_kunderfranco[["gex.logr"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_kunderfranco <- c(mae_kunderfranco, epic = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_kunderfranco, overwrite = TRUE)


# Ren et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_ren[["gex.relz"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_ren <- c(mae_ren, epic = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_ren, overwrite = TRUE)

# Sun et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_sun[["gex.rma"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_sun <- c(curatedPCaData::mae_sun, epic = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_sun, overwrite = TRUE)

# Taylor et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_taylor[["gex.rma"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_taylor <- c(curatedPCaData::mae_taylor, epic = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_taylor, overwrite = TRUE)

# TCGA

tmp <- as.data.frame(immunedeconv::deconvolute(mae_tcga[["gex.tpm"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_tcga <- c(curatedPCaData::mae_tcga, epic = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_tcga, overwrite = TRUE)

# Wang et al.

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_wang[["gex.rma"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_wang <- c(curatedPCaData::mae_wang, epic = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_wang, overwrite = TRUE)

# Kim et al.

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_kim[["gex.rma"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_kim <- c(curatedPCaData::mae_kim, epic = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_kim, overwrite = TRUE)

# Wallace et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_wallace[["gex.rma"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_wallace <- c(mae_wallace, epic = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_wallace, overwrite = TRUE)

# IGC

tmp <- as.data.frame(immunedeconv::deconvolute(mae_igc[["gex.rma"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_igc <- c(mae_igc, epic = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_igc, overwrite = TRUE)

# Weiner et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_weiner[["gex.rma"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_weiner <- c(mae_weiner, epic = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_weiner, overwrite = TRUE)

# True et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_true[["gex.logr"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_true <- c(mae_true, epic = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_true, overwrite = TRUE)


#####################################################
#####################################################
##                                                 ##
##                     QUANTISEQ                   ##
##                                                 ##
#####################################################
#####################################################

# Abida et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_abida[["gex.relz"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_abida <- c(mae_abida, quantiseq = tmp)

# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_abida, overwrite = TRUE)

# Barbieri et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_barbieri[["gex.relz"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_barbieri <- c(mae_barbieri, quantiseq = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_barbieri, overwrite = TRUE)

# Barwick et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_barwick[["gex.logq"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_barwick <- c(mae_barwick, quantiseq = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_barwick, overwrite = TRUE)

# Chandran et al.  


tmp <- as.data.frame(immunedeconv::deconvolute(mae_chandran[["gex.rma"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_chandran <- c(mae_chandran, quantiseq = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_chandran, overwrite = TRUE)

# Friedrich et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_friedrich[["gex.logq"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_friedrich <- c(mae_friedrich, quantiseq = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_friedrich, overwrite = TRUE)

# ICGCCA

tmp <- as.data.frame(immunedeconv::deconvolute(mae_icgcca[["gex.rma"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_icgcca <- c(mae_icgcca, quantiseq = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_icgcca, overwrite = TRUE)

# Kunderfranco et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_kunderfranco[["gex.logr"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_kunderfranco <- c(mae_kunderfranco, quantiseq = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_kunderfranco, overwrite = TRUE)

# Ren et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_ren[["gex.relz"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_ren <- c(mae_ren, quantiseq = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_ren, overwrite = TRUE)

# Sun et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_sun[["gex.rma"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_sun <- c(curatedPCaData::mae_sun, quantiseq = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_sun, overwrite = TRUE)

# Taylor et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_taylor[["gex.rma"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_taylor <- c(curatedPCaData::mae_taylor, quantiseq = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_taylor, overwrite = TRUE)

# TCGA

tmp <- as.data.frame(immunedeconv::deconvolute(mae_tcga[["gex.tpm"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_tcga <- c(curatedPCaData::mae_tcga, quantiseq = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_tcga, overwrite = TRUE)

# Wang et al.

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_wang[["gex.rma"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_wang <- c(curatedPCaData::mae_wang, quantiseq = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_wang, overwrite = TRUE)

# Kim et al.

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_kim[["gex.rma"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_kim <- c(curatedPCaData::mae_kim, quantiseq = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_kim, overwrite = TRUE)

# Wallace et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_wallace[["gex.rma"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_wallace <- c(mae_wallace, quantiseq = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_wallace, overwrite = TRUE)

# IGC

tmp <- as.data.frame(immunedeconv::deconvolute(mae_igc[["gex.rma"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_igc <- c(mae_igc, quantiseq = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_igc, overwrite = TRUE)

# Weiner et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_weiner[["gex.rma"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_weiner <- c(mae_weiner, quantiseq = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_weiner, overwrite = TRUE)

# True et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_true[["gex.logr"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_true <- c(mae_true, quantiseq = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_true, overwrite = TRUE)



#####################################################
#####################################################
##                                                 ##
##                     MCP COUNTER                 ##
##                                                 ##
#####################################################
#####################################################

# Abida et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_abida[["gex.relz"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_abida <- c(mae_abida, mcp = tmp)

# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_abida, overwrite = TRUE)

# Barbieri et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_barbieri[["gex.relz"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_barbieri <- c(mae_barbieri, mcp = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_barbieri, overwrite = TRUE)

# Barwick et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_barwick[["gex.logq"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_barwick <- c(mae_barwick, mcp = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_barwick, overwrite = TRUE)

# Chandran et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_chandran[["gex.rma"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_chandran <- c(mae_chandran, mcp = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_chandran, overwrite = TRUE)

# Friedrich et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_friedrich[["gex.logq"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_friedrich <- c(mae_friedrich, mcp = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_friedrich, overwrite = TRUE)

# ICGCCA

tmp <- as.data.frame(immunedeconv::deconvolute(mae_icgcca[["gex.rma"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_icgcca <- c(mae_icgcca, mcp = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_icgcca, overwrite = TRUE)

# Kunderfranco et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_kunderfranco[["gex.logr"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
#mae_wallace <- c(curatedPCaData::mae_wallace, mcp_counter = tmp)
mae_kunderfranco <- c(curatedPCaData::mae_kunderfranco, mcp = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_kunderfranco, overwrite = TRUE)

# Ren et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_ren[["gex.relz"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_ren <- c(mae_ren, mcp = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_ren, overwrite = TRUE)

# Sun et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_sun[["gex.rma"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_sun <- c(curatedPCaData::mae_sun, mcp = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_sun, overwrite = TRUE)

# Taylor et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_taylor[["gex.rma"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_taylor <- c(curatedPCaData::mae_taylor, mcp = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_taylor, overwrite = TRUE)

# TCGA

tmp <- as.data.frame(immunedeconv::deconvolute(mae_tcga[["gex.tpm"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_tcga <- c(curatedPCaData::mae_tcga, mcp = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_tcga, overwrite = TRUE)

# Wang et al.

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_wang[["gex.rma"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_wang <- c(curatedPCaData::mae_wang, mcp = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_wang, overwrite = TRUE)

# Kim et al.

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_kim[["gex.rma"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_kim <- c(curatedPCaData::mae_kim, mcp = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_kim, overwrite = TRUE)

# Wallace et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_wallace[["gex.rma"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_wallace <- c(mae_wallace, mcp = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_wallace, overwrite = TRUE)

# IGC

tmp <- as.data.frame(immunedeconv::deconvolute(mae_igc[["gex.rma"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_igc <- c(mae_igc, mcp = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_igc, overwrite = TRUE)

# Weiner et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_weiner[["gex.rma"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_weiner <- c(mae_weiner, mcp = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_weiner, overwrite = TRUE)

# True et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_true[["gex.logr"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_true <- c(mae_true, mcp = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_true, overwrite = TRUE)



#####
##
## Genomic risk scores: Prolaris, OncotypeDX & Decipher
## AR scores as used by TCGA, originally presented in Hieronymus et al. 2006
##
#####


## Abida

mae_abida <- c(mae_abida,
	scores = rbind(
		Prolaris = curatedPCaData:::genomic_risk(curatedPCaData::mae_abida, object = "gex_polyA", test = "Prolaris"),
		AR_score = curatedPCaData:::genomic_score(curatedPCaData::mae_abida, object = "gex_polyA", test = "AR")
	)
)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_abida, overwrite = TRUE)

## Barbieri

mae_barbieri <- c(mae_barbieri,
	scores = rbind(
		Prolaris = curatedPCaData:::genomic_risk(curatedPCaData::mae_barbieri, object = "gex", test = "Prolaris"),
		AR_score = curatedPCaData:::genomic_score(curatedPCaData::mae_barbieri, object = "gex", test = "AR")
	)
)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_barbieri, overwrite = TRUE)

## Chandran

mae_chandran <- c(mae_chandran,
	scores = rbind(
		Prolaris = curatedPCaData:::genomic_risk(curatedPCaData::mae_chandran, object = "gex", test = "Prolaris"),
		AR_score = curatedPCaData:::genomic_score(curatedPCaData::mae_chandran, object = "gex", test = "AR")
	)
)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_chandran, overwrite = TRUE)

## Friedrich

mae_friedrich <- c(mae_friedrich,
	scores = rbind(
		Prolaris = curatedPCaData:::genomic_risk(curatedPCaData::mae_friedrich, object = "gex", test = "Prolaris"),
		AR_score = curatedPCaData:::genomic_score(curatedPCaData::mae_friedrich, object = "gex", test = "AR")
	)
)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_friedrich, overwrite = TRUE)

## Hieronymus et al without gex, cannot estimate scores

## ICGC-CA

mae_icgcca <- c(mae_icgcca,
	scores = rbind(
		Prolaris = curatedPCaData:::genomic_risk(curatedPCaData::mae_icgcca, object = "gex", test = "Prolaris"),
		AR_score = curatedPCaData:::genomic_score(curatedPCaData::mae_icgcca, object = "gex", test = "AR")
	)
)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_icgcca, overwrite = TRUE)

## IGC

mae_igc <- c(mae_igc,
	scores = rbind(
		Prolaris = curatedPCaData:::genomic_risk(curatedPCaData::mae_igc, object = "gex", test = "Prolaris"),
		AR_score = curatedPCaData:::genomic_score(curatedPCaData::mae_igc, object = "gex", test = "AR")
	)
)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_igc, overwrite = TRUE)

## Kim

mae_kim <- c(mae_kim,
	scores = rbind(
		Prolaris = curatedPCaData:::genomic_risk(curatedPCaData::mae_kim, object = "gex", test = "Prolaris"),
		AR_score = curatedPCaData:::genomic_score(curatedPCaData::mae_kim, object = "gex", test = "AR")
	)
)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_kim, overwrite = TRUE)

## Kunderfranco

mae_kunderfranco <- c(mae_kunderfranco,
	scores = rbind(
		Prolaris = curatedPCaData:::genomic_risk(curatedPCaData::mae_kunderfranco, object = "gex", test = "Prolaris"),
		AR_score = curatedPCaData:::genomic_score(curatedPCaData::mae_kunderfranco, object = "gex", test = "AR")
	)
)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_kunderfranco, overwrite = TRUE)

## Ren

mae_ren <- c(mae_ren,
	scores = rbind(
		Prolaris = curatedPCaData:::genomic_risk(curatedPCaData::mae_ren, object = "gex", test = "Prolaris"),
		AR_score = curatedPCaData:::genomic_score(curatedPCaData::mae_ren, object = "gex", test = "AR")
	)
)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_ren, overwrite = TRUE)

## Sun

mae_sun <- c(mae_sun,
	scores = rbind(
		Prolaris = curatedPCaData:::genomic_risk(curatedPCaData::mae_sun, object = "gex", test = "Prolaris"),
		AR_score = curatedPCaData:::genomic_score(curatedPCaData::mae_sun, object = "gex", test = "AR")
	)
)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_sun, overwrite = TRUE)

## Taylor

mae_taylor <- c(mae_taylor,
	scores = rbind(
		Prolaris = curatedPCaData:::genomic_risk(curatedPCaData::mae_taylor, object = "gex", test = "Prolaris"),
		AR_score = curatedPCaData:::genomic_score(curatedPCaData::mae_taylor, object = "gex", test = "AR")
	)
)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_taylor, overwrite = TRUE)

## TCGA 

mae_tcga <- c(mae_tcga,
	scores = rbind(
		Prolaris = curatedPCaData:::genomic_risk(curatedPCaData::mae_tcga, object = "gex", test = "Prolaris"),
		AR_score = curatedPCaData:::genomic_score(curatedPCaData::mae_tcga, object = "gex", test = "AR")
	)
)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_tcga, overwrite = TRUE,internal=TRUE)

## True

mae_true <- c(mae_true,
	scores = rbind(
		Prolaris = curatedPCaData:::genomic_risk(curatedPCaData::mae_true, object = "gex", test = "Prolaris"),
		AR_score = curatedPCaData:::genomic_score(curatedPCaData::mae_true, object = "gex", test = "AR")
	)
)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_true, overwrite = TRUE)

## Wallace

mae_wallace <- c(mae_wallace,
	scores = rbind(
		Prolaris = curatedPCaData:::genomic_risk(curatedPCaData::mae_wallace, object = "gex", test = "Prolaris"),
		AR_score = curatedPCaData:::genomic_score(curatedPCaData::mae_wallace, object = "gex", test = "AR")
	)
)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_wallace, overwrite = TRUE)

## IGC

mae_igc <- c(mae_igc,
                 scores = rbind(
                   Prolaris = curatedPCaData:::genomic_risk(curatedPCaData::mae_igc, object = "gex", test = "Prolaris"),
                   AR_score = curatedPCaData:::genomic_score(curatedPCaData::mae_igc, object = "gex", test = "AR")
                 )
)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_igc, overwrite = TRUE)

## Wang

mae_wang <- c(mae_wang,
	scores = rbind(
		Prolaris = curatedPCaData:::genomic_risk(curatedPCaData::mae_wang, object = "gex", test = "Prolaris"),
		AR_score = curatedPCaData:::genomic_score(curatedPCaData::mae_wang, object = "gex", test = "AR")
	)
)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_wang, overwrite = TRUE)

## Weiner

mae_weiner <- c(mae_weiner,
	scores = rbind(
		Prolaris = curatedPCaData:::genomic_risk(curatedPCaData::mae_weiner, object = "gex", test = "Prolaris"),
		AR_score = curatedPCaData:::genomic_score(curatedPCaData::mae_weiner, object = "gex", test = "AR")
	)
)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_weiner, overwrite = TRUE)


####
##
## Tumor purity estimation - DeMix(T) (DeMix is deprecated, DeMixT is the up-to-date version)
##
###

# devtools::install_github("wwylab/DeMixT")

library(DeMixT)


####
##
## Tumor purity estimation - ABSOLUTE
##
###

# Unclear as to how it's available, i.e.

# https://software.broadinstitute.org/cancer/cga/absolute
# (lacking download in https://software.broadinstitute.org/cancer/cga/absolute_run )
# https://software.broadinstitute.org/cancer/cga/absolute_download
# -> Requires separate registration & license agreement, then 'ABSOLUTE_1.0.6.tar.gz' is available

# install.packages("ABSOLUTE_1.0.6.tar.gz", repo=NULL)

# or possibly
# https://www.genepattern.org/modules/docs/ABSOLUTE

library(ABSOLUTE)

# Read in TCGA segmentation file required as ABSOLUTE input
# downloaded from https://www.cbioportal.org/study/cnSegments?id=prad_tcga_pub > "Download a copy number segment file for the selected samples"
TCGA_seg <- read.table("..\\temp\\prad_tcga_pub_segments.seg", sep="\t", header=TRUE)
# ABSOLUTE requires pre-specified column names "Chromosome", "Start", "End", "Num_Probes", and "Segment_Mean";
# replacing default column names with these

