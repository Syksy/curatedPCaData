#####
##
## -- DERIVING NEW VARIABLES --
## Set of scripts that derive new variables from the original raw data
##
#####

## IMMUNE DECONVOLUTION VARIABLES

library(immunedeconv) # NECESSARY!

## Assistance function for adding a new slot in case it doesn't exist, or to avoid existing slot name overlap
addSlotMAE <- function(mae, ...){
	# Transform to a list
	vals <- list(...)
	nams <- names(vals)
	# Iterate over ellipsis slot additions
	for(nam in nams){
		# New slot
		if(!nam %in% names(mae)){
			eval(parse(text=paste("mae <- c(mae,", nam,"=vals[[nam]])", collapse="")))
		# Replace existing slot
		}else{
			mae[[nam]] <- vals[[nam]]
		}
	}
	# Return the MultiAssayExperiment object
	mae
}

## AGGREGATE LOAD
## Load temporary MAE-objects so the whole script can be run and output saved over all iterated 
## Use latest MAE-objects from the package, coming with GEX and other slots necessary for calculating derived variables here-in
##
# Processed datasets
# - Abida et al.
mae_abida <- curatedPCaData::mae_abida
# - Barbieri et al.
mae_barbieri <- curatedPCaData::mae_barbieri
# - Barwick et al.
mae_barwick <- curatedPCaData::mae_barwick
# - Chandran et al.
mae_chandran <- curatedPCaData::mae_chandran
# - Friedrich et al.
mae_friedrich <- curatedPCaData::mae_friedrich
# - ICGC CA
mae_icgcca <- curatedPCaData::mae_icgcca
# - IGC
mae_igc <- curatedPCaData::mae_igc
# - Kim et al.
mae_kim <- curatedPCaData::mae_kim
# - Kunderfranco et al.
mae_kunderfranco <- curatedPCaData::mae_kunderfranco
# - Ren et al.
mae_ren <- curatedPCaData::mae_ren
# - Sun et al.
mae_sun <- curatedPCaData::mae_sun
# - Taylor et al.
mae_taylor <- curatedPCaData::mae_taylor
# - TCGA
mae_tcga <- curatedPCaData::mae_tcga
# - True et al.
mae_true <- curatedPCaData::mae_true
# - Wallace et al.
mae_wallace <- curatedPCaData::mae_wallace
# - Wang et al.
mae_wang <- curatedPCaData::mae_wang
# - Weiner et al.
mae_weiner <- curatedPCaData::mae_weiner

## IMMUNE DECONVOLUTION VARIABLES

#####################################################
#####################################################
##                                                 ##
##                     CIBERSORTX                  ##
##                                                 ##
#####################################################
#####################################################

# Abida et al.
cibersort_polyA_abida<-rio::import("data-raw/CIBERSORTx_abida_Results.csv")
cibersort_polyA_abida<-cibersort_polyA_abida[ , -which(names(cibersort_polyA_abida) %in% c("P-value","Correlation","RMSE"))]
cibersort_polyA_abida <- t(cibersort_polyA_abida)
colnames(cibersort_polyA_abida)<- cibersort_polyA_abida[1,]
cibersort_polyA_abida<-cibersort_polyA_abida[-1,]
cibersort_polyA_abida<-as.matrix(cibersort_polyA_abida)
class(cibersort_polyA_abida) <- "numeric"
#save(cibersort_polyA_abida, file="data-raw/cibersort_polyA_abida.RData")
mae_abida <- addSlotMAE(mae_abida, cibersort = cibersort_polyA_abida)
#mae_capture_abida <- create_mae(study_name = "capture_abida")
#usethis::use_data(mae_abida, overwrite = TRUE)

# Barbieri et al.
cibersort_barbieri<-rio::import("data-raw/CIBERSORTx_barbieri_Results.csv")
cibersort_barbieri<-cibersort_barbieri[ , -which(names(cibersort_barbieri) %in% c("P-value","Correlation","RMSE"))]
cibersort_barbieri <- t(cibersort_barbieri)
colnames(cibersort_barbieri)<- cibersort_barbieri[1,]
cibersort_barbieri<-cibersort_barbieri[-1,]
cibersort_barbieri<-as.matrix(cibersort_barbieri)
class(cibersort_barbieri) <- "numeric"
#save(cibersort_barbieri, file="data-raw/cibersort_barbieri.RData")
mae_barbieri <- addSlotMAE(mae_barbieri, cibersort = cibersort_barbieri)
#mae_barbieri <- create_mae(study_name = "barbieri")
#usethis::use_data(mae_barbieri, overwrite = TRUE)

# Barwick et al.
cibersort_barwick<-rio::import("data-raw/CIBERSORTx_barwick_Results.csv")
cibersort_barwick<-cibersort_barwick[ , -which(names(cibersort_barwick) %in% c("P-value","Correlation","RMSE"))]
cibersort_barwick <- t(cibersort_barwick)
colnames(cibersort_barwick)<- cibersort_barwick[1,]
cibersort_barwick<-cibersort_barwick[-1,]
cibersort_barwick<-as.matrix(cibersort_barwick)
class(cibersort_barwick) <- "numeric"
#save(cibersort_barwick, file="data-raw/cibersort_barwick.RData")
mae_barwick <- addSlotMAE(mae_barwick, cibersort = cibersort_barwick)
#mae_barwick <- create_mae(study_name = "barwick")
#usethis::use_data(mae_barwick, overwrite = TRUE)

# Chandran et al.
cibersort_chandran<-rio::import("data-raw/CIBERSORTx_chandran_Results.csv")
cibersort_chandran<-cibersort_chandran[ , -which(names(cibersort_chandran) %in% c("P-value","Correlation","RMSE"))]
cibersort_chandran <- t(cibersort_chandran)
colnames(cibersort_chandran)<- cibersort_chandran[1,]
cibersort_chandran<-cibersort_chandran[-1,]
cibersort_chandran<-as.matrix(cibersort_chandran)
class(cibersort_chandran) <- "numeric"
#save(cibersort_chandran, file="data-raw/cibersort_chandran.RData")
mae_chandran <- addSlotMAE(mae_chandran, cibersort = cibersort_chandran)
#mae_chandran <- create_mae(study_name = "chandran")
#usethis::use_data(mae_chandran, overwrite = TRUE)

# Friedrich et al.
cibersort_friedrich<-rio::import("data-raw/CIBERSORTx_friedrich_Results.csv")
cibersort_friedrich<-cibersort_friedrich[ , -which(names(cibersort_friedrich) %in% c("P-value","Correlation","RMSE"))]
cibersort_friedrich <- t(cibersort_friedrich)
colnames(cibersort_friedrich)<- cibersort_friedrich[1,]
cibersort_friedrich<-cibersort_friedrich[-1,]
cibersort_friedrich<-as.matrix(cibersort_friedrich)
class(cibersort_friedrich) <- "numeric"
#save(cibersort_friedrich, file="data-raw/cibersort_friedrich.RData")
mae_friedrich <- addSlotMAE(mae_friedrich, cibersort = cibersort_friedrich)
#mae_friedrich <- create_mae(study_name = "friedrich")
#usethis::use_data(mae_friedrich, overwrite = TRUE)

# ICGC Canadian set
cibersort_icgcca<-rio::import("data-raw/CIBERSORTx_icgcca_Results.csv")
cibersort_icgcca<-cibersort_icgcca[ , -which(names(cibersort_icgcca) %in% c("P-value","Correlation","RMSE"))]
cibersort_icgcca <- t(cibersort_icgcca)
colnames(cibersort_icgcca)<- cibersort_icgcca[1,]
cibersort_icgcca<-cibersort_icgcca[-1,]
cibersort_icgcca<-as.matrix(cibersort_icgcca)
class(cibersort_icgcca) <- "numeric"
#save(cibersort_icgcca, file="data-raw/cibersort_icgcca.RData")
#mae_icgcca[["cibersort"]]<-NULL
mae_icgcca <- addSlotMAE(mae_icgcca, cibersort = cibersort_icgcca)
#mae_icgcca <- create_mae(study_name = "icgcca")
#usethis::use_data(mae_icgcca, overwrite = TRUE)

# IGC
cibersort_igc<-rio::import("data-raw/CIBERSORTx_igc_Results.csv")
cibersort_igc<-cibersort_igc[ , -which(names(cibersort_igc) %in% c("P-value","Correlation","RMSE"))]
cibersort_igc <- t(cibersort_igc)
colnames(cibersort_igc)<- cibersort_igc[1,]
cibersort_igc<-cibersort_igc[-1,]
cibersort_igc<-as.matrix(cibersort_igc)
class(cibersort_igc) <- "numeric"
#save(cibersort_igc, file="data-raw/cibersort_igc.RData")
mae_igc <- addSlotMAE(mae_igc, cibersort = cibersort_igc)
#mae_igc <- create_mae(study_name = "igc")
#usethis::use_data(mae_igc, overwrite = TRUE)

# Kim et al.
cibersort_kim<-rio::import("data-raw/CIBERSORTx_kim_Results.csv")
cibersort_kim<-cibersort_kim[ , -which(names(cibersort_kim) %in% c("P-value","Correlation","RMSE"))]
cibersort_kim <- t(cibersort_kim)
colnames(cibersort_kim)<- cibersort_kim[1,]
cibersort_kim<-cibersort_kim[-1,]
cibersort_kim<-as.matrix(cibersort_kim)
class(cibersort_kim) <- "numeric"
#save(cibersort_kim, file="data-raw/cibersort_kim.RData")
mae_kim <- addSlotMAE(mae_kim, cibersort = cibersort_kim)
#mae_kim <- create_mae(study_name = "kim")
#usethis::use_data(mae_kim, overwrite = TRUE)

# Kunderfranco et al.
cibersort_kunderfranco<-rio::import("data-raw/CIBERSORTx_kunderfranco_Results.csv")
cibersort_kunderfranco<-cibersort_kunderfranco[ , -which(names(cibersort_kunderfranco) %in% c("P-value","Correlation","RMSE"))]
cibersort_kunderfranco <- t(cibersort_kunderfranco)
colnames(cibersort_kunderfranco)<- cibersort_kunderfranco[1,]
cibersort_kunderfranco<-cibersort_kunderfranco[-1,]
cibersort_kunderfranco<-as.matrix(cibersort_kunderfranco)
class(cibersort_kunderfranco) <- "numeric"
#save(cibersort_kunderfranco, file="data-raw/cibersort_kunderfranco.RData")
mae_kunderfranco <- addSlotMAE(mae_kunderfranco, cibersort = cibersort_kunderfranco)
#mae_kunderfranco <- create_mae(study_name = "kunderfranco")
#usethis::use_data(mae_kunderfranco, overwrite = TRUE)

# Ren et al.
cibersort_ren<-rio::import("data-raw/CIBERSORTx_ren_Results.csv")
cibersort_ren<-cibersort_ren[ , -which(names(cibersort_ren) %in% c("P-value","Correlation","RMSE"))]
cibersort_ren <- t(cibersort_ren)
colnames(cibersort_ren)<- cibersort_ren[1,]
cibersort_ren<-cibersort_ren[-1,]
cibersort_ren<-as.matrix(cibersort_ren)
class(cibersort_ren) <- "numeric"
#save(cibersort_ren, file="data-raw/cibersort_ren.RData")
mae_ren <- addSlotMAE(mae_ren, cibersort = cibersort_ren)
#mae_ren <- create_mae(study_name = "ren")
#usethis::use_data(mae_ren, overwrite = TRUE)

# Sun et al.
cibersort_sun<-rio::import("data-raw/CIBERSORTx_sun_Results.csv")
cibersort_sun<-cibersort_sun[ , -which(names(cibersort_sun) %in% c("P-value","Correlation","RMSE"))]
cibersort_sun <- t(cibersort_sun)
colnames(cibersort_sun)<- cibersort_sun[1,]
cibersort_sun<-cibersort_sun[-1,]
cibersort_sun<-as.matrix(cibersort_sun)
class(cibersort_sun) <- "numeric"
#save(cibersort_sun, file="data-raw/cibersort_sun.RData")
mae_sun <- addSlotMAE(mae_sun, cibersort = cibersort_sun)
#mae_sun <- create_mae(study_name = "sun")
#usethis::use_data(mae_sun, overwrite = TRUE)

# Taylor et al.
cibersort_taylor<-rio::import("data-raw/CIBERSORTx_taylor_Results.csv")
cibersort_taylor<-cibersort_taylor[ , -which(names(cibersort_taylor) %in% c("P-value","Correlation","RMSE"))]
cibersort_taylor <- t(cibersort_taylor)
colnames(cibersort_taylor)<- cibersort_taylor[1,]
cibersort_taylor<-cibersort_taylor[-1,]
cibersort_taylor<-as.matrix(cibersort_taylor)
class(cibersort_taylor) <- "numeric"
#save(cibersort_taylor, file="data-raw/cibersort_taylor.RData")
mae_taylor <- addSlotMAE(mae_taylor, cibersort = cibersort_taylor)
#mae_taylor <- create_mae(study_name = "taylor")
#usethis::use_data(mae_taylor, overwrite = TRUE)

# TCGA et al.
cibersort_tcga<-rio::import("data-raw/CIBERSORTx_tcga_Results.csv")
cibersort_tcga<-cibersort_tcga[ , -which(names(cibersort_tcga) %in% c("P-value","Correlation","RMSE"))]
cibersort_tcga <- t(cibersort_tcga)
colnames(cibersort_tcga)<- cibersort_tcga[1,]
cibersort_tcga<-cibersort_tcga[-1,]
cibersort_tcga<-as.matrix(cibersort_tcga)
class(cibersort_tcga) <- "numeric"
#save(cibersort_tcga, file="data-raw/cibersort_tcga.RData")
mae_tcga <- addSlotMAE(mae_tcga, cibersort = cibersort_tcga)
#mae_tcga <- create_mae(study_name = "tcga")
#usethis::use_data(mae_tcga, overwrite = TRUE)

# True et al.
cibersort_true<-rio::import("data-raw/CIBERSORTx_true_Results.csv")
cibersort_true<-cibersort_true[ , -which(names(cibersort_true) %in% c("P-value","Correlation","RMSE"))]
cibersort_true <- t(cibersort_true)
colnames(cibersort_true)<- cibersort_true[1,]
cibersort_true<-cibersort_true[-1,]
cibersort_true<-as.matrix(cibersort_true)
class(cibersort_true) <- "numeric"
#save(cibersort_true, file="data-raw/cibersort_true.RData")
mae_true <- addSlotMAE(mae_true, cibersort = cibersort_true)
#mae_true <- create_mae(study_name = "true")
#usethis::use_data(mae_true, overwrite = TRUE)

# Wallace et al.
cibersort_wallace<-rio::import("data-raw/CIBERSORTx_wallace_Results.csv")
cibersort_wallace<-cibersort_wallace[ , -which(names(cibersort_wallace) %in% c("P-value","Correlation","RMSE"))]
cibersort_wallace <- t(cibersort_wallace)
colnames(cibersort_wallace)<- cibersort_wallace[1,]
cibersort_wallace<-cibersort_wallace[-1,]
cibersort_wallace<-as.matrix(cibersort_wallace)
class(cibersort_wallace) <- "numeric"
#save(cibersort_wallace, file="data-raw/cibersort_wallace.RData")
mae_wallace <- addSlotMAE(mae_wallace, cibersort = cibersort_wallace)
#mae_wallace <- create_mae(study_name = "wallace")
#usethis::use_data(mae_wallace, overwrite = TRUE)

# Wang et al.
cibersort_wang<-rio::import("data-raw/CIBERSORTx_wang_Results.csv")
cibersort_wang<-cibersort_wang[ , -which(names(cibersort_wang) %in% c("P-value","Correlation","RMSE"))]
cibersort_wang <- t(cibersort_wang)
colnames(cibersort_wang)<- cibersort_wang[1,]
cibersort_wang<-cibersort_wang[-1,]
cibersort_wang<-as.matrix(cibersort_wang)
class(cibersort_wang) <- "numeric"
#save(cibersort_wang, file="data-raw/cibersort_wang.RData")
mae_wang <- addSlotMAE(mae_wang, cibersort = cibersort_wang)
#mae_wang <- create_mae(study_name = "wang")
#usethis::use_data(mae_wang, overwrite = TRUE)

# Weiner et al.
cibersort_weiner<-rio::import("data-raw/CIBERSORTx_weiner_Results.csv")
cibersort_weiner<-cibersort_weiner[ , -which(names(cibersort_weiner) %in% c("P-value","Pearson Correlation","RMSE","Absolute score"))]
cibersort_weiner <- t(cibersort_weiner)
colnames(cibersort_weiner)<- cibersort_weiner[1,]
cibersort_weiner<-cibersort_weiner[-1,]
cibersort_weiner<-as.matrix(cibersort_weiner)
class(cibersort_weiner) <- "numeric"
#save(cibersort_weiner, file="data-raw/cibersort_weiner.RData")
mae_weiner <- addSlotMAE(mae_weiner, cibersort = cibersort_weiner)

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
mae_abida <- addSlotMAE(mae_abida, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_abida, overwrite = TRUE)

# Barbieri et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_barbieri[["gex.relz"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_barbieri <- addSlotMAE(mae_barbieri, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_barbieri, overwrite = TRUE)

## METHOD FAILS DUE TO LACK OF GENE OVERLAP
if(FALSE){
	# Barwick et al.
	tmp <- as.data.frame(immunedeconv::deconvolute(mae_barwick[["gex.logq"]], method="xcell"))
	rownames(tmp) <- tmp$cell_type
	# Omit cell type column and store only data of cell type populations
	tmp <- as.matrix(tmp[,-1])
	mae_barwick <- addSlotMAE(mae_barwick, xcell = tmp)
	# Save the derived new 'assay' types to the mae-object
	#usethis::use_data(mae_barwick, overwrite = TRUE)
}

# Chandran et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_chandran[["gex.rma"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_chandran <- addSlotMAE(mae_chandran, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_chandran, overwrite = TRUE)

# Friedrich et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_friedrich[["gex.logq"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_friedrich <- addSlotMAE(mae_friedrich, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_friedrich, overwrite = TRUE)

# ICGCCA
tmp <- as.data.frame(immunedeconv::deconvolute(mae_icgcca[["gex.rma"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_icgcca <- addSlotMAE(mae_icgcca, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_icgcca, overwrite = TRUE)

# Kunderfranco et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_kunderfranco[["gex.logr"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_kunderfranco <- addSlotMAE(mae_kunderfranco, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_kunderfranco, overwrite = TRUE)

# Ren et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_ren[["gex.relz"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_ren <- addSlotMAE(mae_ren, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_ren, overwrite = TRUE)

# Sun et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_sun[["gex.rma"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_sun <- addSlotMAE(mae_sun, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_sun, overwrite = TRUE)

# Taylor et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_taylor[["gex.rma"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_taylor <- addSlotMAE(curatedPCaData::mae_taylor, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_taylor, overwrite = TRUE)

# TCGA
tmp <- as.data.frame(immunedeconv::deconvolute(mae_tcga[["gex.fpkm"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_tcga <- addSlotMAE(mae_tcga, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_tcga, overwrite = TRUE)

# Wang et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_wang[["gex.rma"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_wang <- addSlotMAE(mae_wang, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_wang, overwrite = TRUE)

# Kim et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_kim[["gex.rma"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_kim <- addSlotMAE(mae_kim, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_kim, overwrite = TRUE)

# Wallace et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_wallace[["gex.rma"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_wallace <- addSlotMAE(mae_wallace, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_wallace, overwrite = TRUE)

# Igc 
tmp <- as.data.frame(immunedeconv::deconvolute(mae_igc[["gex.rma"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_igc <- addSlotMAE(mae_igc, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_igc, overwrite = TRUE)

# Weiner et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_weiner[["gex.rma"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_weiner <- addSlotMAE(mae_weiner, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_weiner, overwrite = TRUE)

## METHOD FAILS DUE TO LACK OF GENE OVERLAP
if(FALSE){
	# True et al.
	tmp <- as.data.frame(immunedeconv::deconvolute(mae_true[["gex.logr"]], method="xcell"))
	rownames(tmp) <- tmp$cell_type
	# Omit cell type column and store only data of cell type populations
	tmp <- as.matrix(tmp[,-1])
	mae_true <- addSlotMAE(mae_true, xcell = tmp)
	# Save the derived new 'assay' types to the mae-object
	#usethis::use_data(mae_true, overwrite = TRUE)
}


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
mae_abida <- addSlotMAE(mae_abida, epic = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_abida, overwrite = TRUE)

# FAILS DUE TO OVERLAP OF GENES
if(FALSE){
	# Barwick et al.
	tmp <- as.data.frame(immunedeconv::deconvolute(mae_barwick[["gex.logq"]], method="epic"))
	rownames(tmp) <- tmp$cell_type
	# Omit cell type column and store only data of cell type populations
	tmp <- as.matrix(tmp[,-1])
	mae_barwick <- addSlotMAE(mae_barwick, epic = tmp)
	# Save the derived new 'assay' types to the mae-object
	#usethis::use_data(mae_barwick, overwrite = TRUE)
}

# Chandran et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_chandran[["gex.rma"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_chandran <- addSlotMAE(mae_chandran, epic = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_chandran, overwrite = TRUE)

# Friedrich et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_friedrich[["gex.logq"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_friedrich <- addSlotMAE(mae_friedrich, epic = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_friedrich, overwrite = TRUE)

# ICGCCA
tmp <- as.data.frame(immunedeconv::deconvolute(mae_icgcca[["gex.rma"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_icgcca <- addSlotMAE(mae_icgcca, epic = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_icgcca, overwrite = TRUE)

# Kunderfranco et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_kunderfranco[["gex.logr"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_kunderfranco <- addSlotMAE(mae_kunderfranco, epic = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_kunderfranco, overwrite = TRUE)

# Ren et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_ren[["gex.relz"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_ren <- addSlotMAE(mae_ren, epic = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_ren, overwrite = TRUE)

# Sun et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_sun[["gex.rma"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_sun <- addSlotMAE(mae_sun, epic = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_sun, overwrite = TRUE)

# Taylor et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_taylor[["gex.rma"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_taylor <- addSlotMAE(mae_taylor, epic = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_taylor, overwrite = TRUE)

# TCGA
tmp <- as.data.frame(immunedeconv::deconvolute(mae_tcga[["gex.fpkm"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_tcga <- addSlotMAE(mae_tcga, epic = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_tcga, overwrite = TRUE)

# Wang et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_wang[["gex.rma"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_wang <- addSlotMAE(mae_wang, epic = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_wang, overwrite = TRUE)

# Kim et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_kim[["gex.rma"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_kim <- addSlotMAE(mae_kim, epic = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_kim, overwrite = TRUE)

# Wallace et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_wallace[["gex.rma"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_wallace <- addSlotMAE(mae_wallace, epic = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_wallace, overwrite = TRUE)

# IGC
tmp <- as.data.frame(immunedeconv::deconvolute(mae_igc[["gex.rma"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_igc <- addSlotMAE(mae_igc, epic = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_igc, overwrite = TRUE)

# Wang et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_wang[["gex.rma"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_wang <- addSlotMAE(mae_wang, epic = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_weiner, overwrite = TRUE)

# Weiner et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_weiner[["gex.rma"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_weiner <- addSlotMAE(mae_weiner, epic = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_weiner, overwrite = TRUE)

# True et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_true[["gex.logr"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_true <- addSlotMAE(mae_true, epic = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_true, overwrite = TRUE)


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
mae_abida <- addSlotMAE(mae_abida, quantiseq = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_abida, overwrite = TRUE)

# Barbieri et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_barbieri[["gex.relz"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_barbieri <- addSlotMAE(mae_barbieri, quantiseq = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_barbieri, overwrite = TRUE)

# Barwick et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_barwick[["gex.logq"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_barwick <- addSlotMAE(mae_barwick, quantiseq = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_barwick, overwrite = TRUE)

# Chandran et al.  
tmp <- as.data.frame(immunedeconv::deconvolute(mae_chandran[["gex.rma"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_chandran <- addSlotMAE(mae_chandran, quantiseq = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_chandran, overwrite = TRUE)

# Friedrich et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_friedrich[["gex.logq"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_friedrich <- addSlotMAE(mae_friedrich, quantiseq = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_friedrich, overwrite = TRUE)

# ICGCCA
tmp <- as.data.frame(immunedeconv::deconvolute(mae_icgcca[["gex.rma"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_icgcca <- addSlotMAE(mae_icgcca, quantiseq = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_icgcca, overwrite = TRUE)

# Kunderfranco et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_kunderfranco[["gex.logr"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_kunderfranco <- addSlotMAE(mae_kunderfranco, quantiseq = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_kunderfranco, overwrite = TRUE)

# Ren et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_ren[["gex.relz"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_ren <- addSlotMAE(mae_ren, quantiseq = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_ren, overwrite = TRUE)

# Sun et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_sun[["gex.rma"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_sun <- addSlotMAE(mae_sun, quantiseq = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_sun, overwrite = TRUE)

# Taylor et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_taylor[["gex.rma"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_taylor <- addSlotMAE(mae_taylor, quantiseq = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_taylor, overwrite = TRUE)

# TCGA
tmp <- as.data.frame(immunedeconv::deconvolute(mae_tcga[["gex.fpkm"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_tcga <- addSlotMAE(mae_tcga, quantiseq = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_tcga, overwrite = TRUE)

# Wang et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_wang[["gex.rma"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_wang <- addSlotMAE(mae_wang, quantiseq = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_wang, overwrite = TRUE)

# Kim et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_kim[["gex.rma"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_kim <- addSlotMAE(mae_kim, quantiseq = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_kim, overwrite = TRUE)

# Wallace et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_wallace[["gex.rma"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_wallace <- addSlotMAE(mae_wallace, quantiseq = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_wallace, overwrite = TRUE)

# IGC
tmp <- as.data.frame(immunedeconv::deconvolute(mae_igc[["gex.rma"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_igc <- addSlotMAE(mae_igc, quantiseq = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_igc, overwrite = TRUE)

# Wang et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_wang[["gex.rma"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_wang <- addSlotMAE(mae_wang, quantiseq = tmp)

# Weiner et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_weiner[["gex.rma"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_weiner <- addSlotMAE(mae_weiner, quantiseq = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_weiner, overwrite = TRUE)

# True et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_true[["gex.logr"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_true <- addSlotMAE(mae_true, quantiseq = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_true, overwrite = TRUE)



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
mae_abida <- addSlotMAE(mae_abida, mcp = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_abida, overwrite = TRUE)

# Barbieri et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_barbieri[["gex.relz"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_barbieri <- addSlotMAE(mae_barbieri, mcp = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_barbieri, overwrite = TRUE)

# Barwick et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_barwick[["gex.logq"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_barwick <- addSlotMAE(mae_barwick, mcp = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_barwick, overwrite = TRUE)

# Chandran et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_chandran[["gex.rma"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_chandran <- addSlotMAE(mae_chandran, mcp = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_chandran, overwrite = TRUE)

# Friedrich et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_friedrich[["gex.logq"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_friedrich <- addSlotMAE(mae_friedrich, mcp = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_friedrich, overwrite = TRUE)

# ICGCCA
tmp <- as.data.frame(immunedeconv::deconvolute(mae_icgcca[["gex.rma"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_icgcca <- addSlotMAE(mae_icgcca, mcp = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_icgcca, overwrite = TRUE)

# Kunderfranco et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_kunderfranco[["gex.logr"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
#mae_wallace <- c(curatedPCaData::mae_wallace, mcp_counter = tmp)
mae_kunderfranco <- addSlotMAE(mae_kunderfranco, mcp = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_kunderfranco, overwrite = TRUE)

# Ren et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_ren[["gex.relz"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_ren <- addSlotMAE(mae_ren, mcp = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_ren, overwrite = TRUE)

# Sun et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_sun[["gex.rma"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_sun <- addSlotMAE(mae_sun, mcp = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_sun, overwrite = TRUE)

# Taylor et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_taylor[["gex.rma"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_taylor <- addSlotMAE(mae_taylor, mcp = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_taylor, overwrite = TRUE)

# TCGA
tmp <- as.data.frame(immunedeconv::deconvolute(mae_tcga[["gex.fpkm"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_tcga <- addSlotMAE(mae_tcga, mcp = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_tcga, overwrite = TRUE)

# Wang et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_wang[["gex.rma"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_wang <- addSlotMAE(mae_wang, mcp = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_wang, overwrite = TRUE)

# Kim et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_kim[["gex.rma"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_kim <- addSlotMAE(mae_kim, mcp = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_kim, overwrite = TRUE)

# Wallace et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_wallace[["gex.rma"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_wallace <- addSlotMAE(mae_wallace, mcp = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_wallace, overwrite = TRUE)

# IGC
tmp <- as.data.frame(immunedeconv::deconvolute(mae_igc[["gex.rma"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_igc <- addSlotMAE(mae_igc, mcp = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_igc, overwrite = TRUE)

# Weiner et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_weiner[["gex.rma"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_weiner <- addSlotMAE(mae_weiner, mcp = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_weiner, overwrite = TRUE)

# True et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_true[["gex.logr"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_true <- addSlotMAE(mae_true, mcp = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_true, overwrite = TRUE)



#####
##
## Genomic risk scores: Prolaris, OncotypeDX & Decipher
## AR scores as used by TCGA, originally presented in Hieronymus et al. 2006
##
#####


## Abida
mae_abida <- addSlotMAE(mae_abida,
	scores = rbind(
		Prolaris = curatedPCaData:::genomic_risk(mae_abida, slot = "gex.relz", test = "Prolaris"),
		AR_score = curatedPCaData:::genomic_score(mae_abida, slot = "gex.relz", test = "AR")
	)
)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_abida, overwrite = TRUE)

## Barbieri
mae_barbieri <- addSlotMAE(mae_barbieri,
	scores = rbind(
		Prolaris = curatedPCaData:::genomic_risk(mae_barbieri, slot = "gex.relz", test = "Prolaris"),
		AR_score = curatedPCaData:::genomic_score(mae_barbieri, slot = "gex.relz", test = "AR")
	)
)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_barbieri, overwrite = TRUE)

## Chandran
mae_chandran <- addSlotMAE(mae_chandran,
	scores = rbind(
		Prolaris = curatedPCaData:::genomic_risk(mae_chandran, slot = "gex.rma", test = "Prolaris"),
		AR_score = curatedPCaData:::genomic_score(mae_chandran, slot = "gex.rma", test = "AR")
	)
)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_chandran, overwrite = TRUE)

## Friedrich
mae_friedrich <- addSlotMAE(mae_friedrich,
	scores = rbind(
		Prolaris = curatedPCaData:::genomic_risk(mae_friedrich, slot = "gex.logq", test = "Prolaris"),
		AR_score = curatedPCaData:::genomic_score(mae_friedrich, slot = "gex.logq", test = "AR")
	)
)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_friedrich, overwrite = TRUE)

## ICGC-CA
mae_icgcca <- addSlotMAE(mae_icgcca,
	scores = rbind(
		Prolaris = curatedPCaData:::genomic_risk(mae_icgcca, slot = "gex.rma", test = "Prolaris"),
		AR_score = curatedPCaData:::genomic_score(mae_icgcca, slot = "gex.rma", test = "AR")
	)
)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_icgcca, overwrite = TRUE)

## IGC
mae_igc <- addSlotMAE(mae_igc,
	scores = rbind(
		Prolaris = curatedPCaData:::genomic_risk(mae_igc, slot = "gex.rma", test = "Prolaris"),
		AR_score = curatedPCaData:::genomic_score(mae_igc, slot = "gex.rma", test = "AR")
	)
)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_igc, overwrite = TRUE)

## Kim
mae_kim <- addSlotMAE(mae_kim,
	scores = rbind(
		Prolaris = curatedPCaData:::genomic_risk(mae_kim, slot = "gex.rma", test = "Prolaris"),
		AR_score = curatedPCaData:::genomic_score(mae_kim, slot = "gex.rma", test = "AR")
	)
)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_kim, overwrite = TRUE)

## Kunderfranco
mae_kunderfranco <- addSlotMAE(mae_kunderfranco,
	scores = rbind(
		Prolaris = curatedPCaData:::genomic_risk(mae_kunderfranco, slot = "gex.logr", test = "Prolaris"),
		AR_score = curatedPCaData:::genomic_score(mae_kunderfranco, slot = "gex.logr", test = "AR")
	)
)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_kunderfranco, overwrite = TRUE)

## Ren
mae_ren <- addSlotMAE(mae_ren,
	scores = rbind(
		Prolaris = curatedPCaData:::genomic_risk(mae_ren, slot = "gex.relz", test = "Prolaris"),
		AR_score = curatedPCaData:::genomic_score(mae_ren, slot = "gex.relz", test = "AR")
	)
)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_ren, overwrite = TRUE)

## Sun
mae_sun <- addSlotMAE(mae_sun,
	scores = rbind(
		Prolaris = curatedPCaData:::genomic_risk(mae_sun, slot = "gex.rma", test = "Prolaris"),
		AR_score = curatedPCaData:::genomic_score(mae_sun, slot = "gex.rma", test = "AR")
	)
)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_sun, overwrite = TRUE)

## Taylor
mae_taylor <- addSlotMAE(mae_taylor,
	scores = rbind(
		Prolaris = curatedPCaData:::genomic_risk(mae_taylor, slot = "gex.rma", test = "Prolaris"),
		AR_score = curatedPCaData:::genomic_score(mae_taylor, slot = "gex.rma", test = "AR")
	)
)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_taylor, overwrite = TRUE)

## TCGA 
mae_tcga <- addSlotMAE(mae_tcga,
	scores = rbind(
		Prolaris = curatedPCaData:::genomic_risk(mae_tcga, slot = "gex.fpkm", test = "Prolaris"),
		AR_score = curatedPCaData:::genomic_score(mae_tcga, slot = "gex.fpkm", test = "AR")
	)
)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_tcga, overwrite = TRUE,internal=TRUE)

## True
mae_true <- addSlotMAE(mae_true,
	scores = rbind(
		Prolaris = curatedPCaData:::genomic_risk(mae_true, slot = "gex.logr", test = "Prolaris"),
		AR_score = curatedPCaData:::genomic_score(mae_true, slot = "gex.logr", test = "AR")
	)
)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_true, overwrite = TRUE)

## Wallace
mae_wallace <- addSlotMAE(mae_wallace,
	scores = rbind(
		Prolaris = curatedPCaData:::genomic_risk(mae_wallace, slot = "gex.rma", test = "Prolaris"),
		AR_score = curatedPCaData:::genomic_score(mae_wallace, slot = "gex.rma", test = "AR")
	)
)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_wallace, overwrite = TRUE)

## IGC
mae_igc <- addSlotMAE(mae_igc,
                 scores = rbind(
                   Prolaris = curatedPCaData:::genomic_risk(mae_igc, slot = "gex.rma", test = "Prolaris"),
                   AR_score = curatedPCaData:::genomic_score(mae_igc, slot = "gex.rma", test = "AR")
                 )
)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_igc, overwrite = TRUE)

## Wang
mae_wang <- addSlotMAE(mae_wang,
	scores = rbind(
		Prolaris = curatedPCaData:::genomic_risk(mae_wang, slot = "gex.rma", test = "Prolaris"),
		AR_score = curatedPCaData:::genomic_score(mae_wang, slot = "gex.rma", test = "AR")
	)
)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_wang, overwrite = TRUE)

## Weiner
mae_weiner <- addSlotMAE(mae_weiner,
	scores = rbind(
		Prolaris = curatedPCaData:::genomic_risk(mae_weiner, slot = "gex.rma", test = "Prolaris"),
		AR_score = curatedPCaData:::genomic_score(mae_weiner, slot = "gex.rma", test = "AR")
	)
)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_weiner, overwrite = TRUE)


## ----

## AGGREGATE SAVE
## One by one save back all the MAE objects back so all updated slots are correctly stored in the new MAE objects
##

# Processed datasets
# - Abida et al.
usethis::use_data(mae_abida, overwrite = TRUE)
# - Barbieri et al.
usethis::use_data(mae_barbieri, overwrite = TRUE)
# - Barwick et al.
usethis::use_data(mae_barwick, overwrite = TRUE)
# - Chandran et al.
usethis::use_data(mae_chandran, overwrite = TRUE)
# - Friedrich et al.
usethis::use_data(mae_friedrich, overwrite = TRUE)
# - ICGC CA
usethis::use_data(mae_icgcca, overwrite = TRUE)
# - IGC
usethis::use_data(mae_igc, overwrite = TRUE)
# - Kim et al.
usethis::use_data(mae_kim, overwrite = TRUE)
# - Kunderfranco et al.
usethis::use_data(mae_kunderfranco, overwrite = TRUE)
# - Ren et al.
usethis::use_data(mae_ren, overwrite = TRUE)
# - Sun et al.
usethis::use_data(mae_sun, overwrite = TRUE)
# - Taylor et al.
usethis::use_data(mae_taylor, overwrite = TRUE)
# - TCGA
usethis::use_data(mae_tcga, overwrite = TRUE)
# - True et al.
usethis::use_data(mae_true, overwrite = TRUE)
# - Wallace et al.
usethis::use_data(mae_wallace, overwrite = TRUE)
# - Wang et al.
usethis::use_data(mae_wang, overwrite = TRUE)
# - Weiner et al.
usethis::use_data(mae_weiner, overwrite = TRUE)

