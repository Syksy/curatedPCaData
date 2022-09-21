#####
##
## -- DERIVING NEW VARIABLES --
## Set of scripts that derive new variables from the original raw data
##
#####

## IMMUNE DECONVOLUTION VARIABLES

library(immunedeconv) # NECESSARY!
library(RaggedExperiment) # Necessary, otherwise MAE objects cannot be fully opened correctly

## Assistance function for adding a new slot in case it doesn't exist, or to avoid existing slot name overlap
addSlotMAE <- function(
	mae,	# MultiAssayExperiment object
	...	# Added slots
){
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
## Update: load the data from working directory '/data/mae_*.rda'
# Processed datasets
# - Abida et al.
load("data/mae_abida.rda")
# - Barbieri et al.
load("data/mae_barbieri.rda")
# - Barwick et al.
load("data/mae_barwick.rda")
# - Chandran et al.
load("data/mae_chandran.rda")
# - Friedrich et al.
load("data/mae_friedrich.rda")
# - ICGC CA
load("data/mae_icgcca.rda")
# - IGC
load("data/mae_igc.rda")
# - Kim et al.
load("data/mae_kim.rda")
# - Kunderfranco et al.
load("data/mae_kunderfranco.rda")
# - Ren et al.
load("data/mae_ren.rda")
# - Sun et al.
load("data/mae_sun.rda")
# - Taylor et al.
load("data/mae_taylor.rda")
# - TCGA
load("data/mae_tcga.rda")
# - True et al.
load("data/mae_true.rda")
# - Wallace et al.
load("data/mae_wallace.rda")
# - Wang et al.
load("data/mae_wang.rda")
# - Weiner et al.
load("data/mae_weiner.rda")

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
mae_abida <- addSlotMAE(mae_abida, cibersort = cibersort_polyA_abida)

# Barbieri et al.
cibersort_barbieri<-rio::import("data-raw/CIBERSORTx_barbieri_Results.csv")
cibersort_barbieri<-cibersort_barbieri[ , -which(names(cibersort_barbieri) %in% c("P-value","Correlation","RMSE"))]
cibersort_barbieri <- t(cibersort_barbieri)
colnames(cibersort_barbieri)<- cibersort_barbieri[1,]
cibersort_barbieri<-cibersort_barbieri[-1,]
cibersort_barbieri<-as.matrix(cibersort_barbieri)
class(cibersort_barbieri) <- "numeric"
mae_barbieri <- addSlotMAE(mae_barbieri, cibersort = cibersort_barbieri)

# Barwick et al.
cibersort_barwick<-rio::import("data-raw/CIBERSORTx_barwick_Results.csv")
cibersort_barwick<-cibersort_barwick[ , -which(names(cibersort_barwick) %in% c("P-value","Correlation","RMSE"))]
cibersort_barwick <- t(cibersort_barwick)
colnames(cibersort_barwick)<- cibersort_barwick[1,]
cibersort_barwick<-cibersort_barwick[-1,]
cibersort_barwick<-as.matrix(cibersort_barwick)
class(cibersort_barwick) <- "numeric"
mae_barwick <- addSlotMAE(mae_barwick, cibersort = cibersort_barwick)

# Chandran et al.
cibersort_chandran<-rio::import("data-raw/CIBERSORTx_chandran_Results.csv")
cibersort_chandran<-cibersort_chandran[ , -which(names(cibersort_chandran) %in% c("P-value","Correlation","RMSE"))]
cibersort_chandran <- t(cibersort_chandran)
colnames(cibersort_chandran)<- cibersort_chandran[1,]
cibersort_chandran<-cibersort_chandran[-1,]
cibersort_chandran<-as.matrix(cibersort_chandran)
class(cibersort_chandran) <- "numeric"
mae_chandran <- addSlotMAE(mae_chandran, cibersort = cibersort_chandran)

# Friedrich et al.
cibersort_friedrich<-rio::import("data-raw/CIBERSORTx_friedrich_Results.csv")
cibersort_friedrich<-cibersort_friedrich[ , -which(names(cibersort_friedrich) %in% c("P-value","Correlation","RMSE"))]
cibersort_friedrich <- t(cibersort_friedrich)
colnames(cibersort_friedrich)<- cibersort_friedrich[1,]
cibersort_friedrich<-cibersort_friedrich[-1,]
cibersort_friedrich<-as.matrix(cibersort_friedrich)
class(cibersort_friedrich) <- "numeric"
mae_friedrich <- addSlotMAE(mae_friedrich, cibersort = cibersort_friedrich)

# ICGC Canadian set
cibersort_icgcca<-rio::import("data-raw/CIBERSORTx_icgcca_Results.csv")
cibersort_icgcca<-cibersort_icgcca[ , -which(names(cibersort_icgcca) %in% c("P-value","Correlation","RMSE"))]
cibersort_icgcca <- t(cibersort_icgcca)
colnames(cibersort_icgcca)<- cibersort_icgcca[1,]
cibersort_icgcca<-cibersort_icgcca[-1,]
cibersort_icgcca<-as.matrix(cibersort_icgcca)
class(cibersort_icgcca) <- "numeric"
mae_icgcca <- addSlotMAE(mae_icgcca, cibersort = cibersort_icgcca)

# IGC
cibersort_igc<-rio::import("data-raw/CIBERSORTx_igc_Results.csv")
cibersort_igc<-cibersort_igc[ , -which(names(cibersort_igc) %in% c("P-value","Correlation","RMSE"))]
cibersort_igc <- t(cibersort_igc)
colnames(cibersort_igc)<- cibersort_igc[1,]
cibersort_igc<-cibersort_igc[-1,]
cibersort_igc<-as.matrix(cibersort_igc)
class(cibersort_igc) <- "numeric"
mae_igc <- addSlotMAE(mae_igc, cibersort = cibersort_igc)

# Kim et al.
cibersort_kim<-rio::import("data-raw/CIBERSORTx_kim_Results.csv")
cibersort_kim<-cibersort_kim[ , -which(names(cibersort_kim) %in% c("P-value","Correlation","RMSE"))]
cibersort_kim <- t(cibersort_kim)
colnames(cibersort_kim)<- cibersort_kim[1,]
cibersort_kim<-cibersort_kim[-1,]
cibersort_kim<-as.matrix(cibersort_kim)
class(cibersort_kim) <- "numeric"
mae_kim <- addSlotMAE(mae_kim, cibersort = cibersort_kim)

# Kunderfranco et al.
cibersort_kunderfranco<-rio::import("data-raw/CIBERSORTx_kunderfranco_Results.csv")
cibersort_kunderfranco<-cibersort_kunderfranco[ , -which(names(cibersort_kunderfranco) %in% c("P-value","Correlation","RMSE"))]
cibersort_kunderfranco <- t(cibersort_kunderfranco)
colnames(cibersort_kunderfranco)<- cibersort_kunderfranco[1,]
cibersort_kunderfranco<-cibersort_kunderfranco[-1,]
cibersort_kunderfranco<-as.matrix(cibersort_kunderfranco)
class(cibersort_kunderfranco) <- "numeric"
mae_kunderfranco <- addSlotMAE(mae_kunderfranco, cibersort = cibersort_kunderfranco)

# Ren et al.
cibersort_ren<-rio::import("data-raw/CIBERSORTx_ren_Results.csv")
cibersort_ren<-cibersort_ren[ , -which(names(cibersort_ren) %in% c("P-value","Correlation","RMSE"))]
cibersort_ren <- t(cibersort_ren)
colnames(cibersort_ren)<- cibersort_ren[1,]
cibersort_ren<-cibersort_ren[-1,]
cibersort_ren<-as.matrix(cibersort_ren)
class(cibersort_ren) <- "numeric"
mae_ren <- addSlotMAE(mae_ren, cibersort = cibersort_ren)

# Sun et al.
cibersort_sun<-rio::import("data-raw/CIBERSORTx_sun_Results.csv")
cibersort_sun<-cibersort_sun[ , -which(names(cibersort_sun) %in% c("P-value","Correlation","RMSE"))]
cibersort_sun <- t(cibersort_sun)
colnames(cibersort_sun)<- cibersort_sun[1,]
cibersort_sun<-cibersort_sun[-1,]
cibersort_sun<-as.matrix(cibersort_sun)
class(cibersort_sun) <- "numeric"
mae_sun <- addSlotMAE(mae_sun, cibersort = cibersort_sun)

# Taylor et al.
cibersort_taylor<-rio::import("data-raw/CIBERSORTx_taylor_Results.csv")
cibersort_taylor<-cibersort_taylor[ , -which(names(cibersort_taylor) %in% c("P-value","Correlation","RMSE"))]
cibersort_taylor <- t(cibersort_taylor)
colnames(cibersort_taylor)<- cibersort_taylor[1,]
cibersort_taylor<-cibersort_taylor[-1,]
cibersort_taylor<-as.matrix(cibersort_taylor)
class(cibersort_taylor) <- "numeric"
mae_taylor <- addSlotMAE(mae_taylor, cibersort = cibersort_taylor)

# TCGA
cibersort_tcga<-rio::import("data-raw/CIBERSORTx_tcga_Results.csv")
cibersort_tcga<-cibersort_tcga[ , -which(names(cibersort_tcga) %in% c("P-value","Correlation","RMSE"))]
cibersort_tcga <- t(cibersort_tcga)
colnames(cibersort_tcga)<- cibersort_tcga[1,]
cibersort_tcga<-cibersort_tcga[-1,]
cibersort_tcga<-as.matrix(cibersort_tcga)
class(cibersort_tcga) <- "numeric"
mae_tcga <- addSlotMAE(mae_tcga, cibersort = cibersort_tcga)

# True et al.
cibersort_true<-rio::import("data-raw/CIBERSORTx_true_Results.csv")
cibersort_true<-cibersort_true[ , -which(names(cibersort_true) %in% c("P-value","Correlation","RMSE"))]
cibersort_true <- t(cibersort_true)
colnames(cibersort_true)<- cibersort_true[1,]
cibersort_true<-cibersort_true[-1,]
cibersort_true<-as.matrix(cibersort_true)
class(cibersort_true) <- "numeric"
mae_true <- addSlotMAE(mae_true, cibersort = cibersort_true)

# Wallace et al.
cibersort_wallace<-rio::import("data-raw/CIBERSORTx_wallace_Results.csv")
cibersort_wallace<-cibersort_wallace[ , -which(names(cibersort_wallace) %in% c("P-value","Correlation","RMSE"))]
cibersort_wallace <- t(cibersort_wallace)
colnames(cibersort_wallace)<- cibersort_wallace[1,]
cibersort_wallace<-cibersort_wallace[-1,]
cibersort_wallace<-as.matrix(cibersort_wallace)
class(cibersort_wallace) <- "numeric"
mae_wallace <- addSlotMAE(mae_wallace, cibersort = cibersort_wallace)

# Wang et al.
cibersort_wang<-rio::import("data-raw/CIBERSORTx_wang_Results.csv")
cibersort_wang<-cibersort_wang[ , -which(names(cibersort_wang) %in% c("P-value","Correlation","RMSE"))]
cibersort_wang <- t(cibersort_wang)
colnames(cibersort_wang)<- cibersort_wang[1,]
cibersort_wang<-cibersort_wang[-1,]
cibersort_wang<-as.matrix(cibersort_wang)
class(cibersort_wang) <- "numeric"
mae_wang <- addSlotMAE(mae_wang, cibersort = cibersort_wang)

# Weiner et al.
cibersort_weiner<-rio::import("data-raw/CIBERSORTx_weiner_Results.csv")
cibersort_weiner<-cibersort_weiner[ , -which(names(cibersort_weiner) %in% c("P-value","Pearson Correlation","RMSE","Absolute score"))]
cibersort_weiner <- t(cibersort_weiner)
colnames(cibersort_weiner)<- cibersort_weiner[1,]
cibersort_weiner<-cibersort_weiner[-1,]
cibersort_weiner<-as.matrix(cibersort_weiner)
class(cibersort_weiner) <- "numeric"
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

# Barbieri et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_barbieri[["gex.relz"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_barbieri <- addSlotMAE(mae_barbieri, xcell = tmp)

## METHOD FAILS DUE TO LACK OF GENE OVERLAP
if(FALSE){
	# Barwick et al.
	tmp <- as.data.frame(immunedeconv::deconvolute(mae_barwick[["gex.logq"]], method="xcell"))
	rownames(tmp) <- tmp$cell_type
	# Omit cell type column and store only data of cell type populations
	tmp <- as.matrix(tmp[,-1])
	mae_barwick <- addSlotMAE(mae_barwick, xcell = tmp)
}

# Chandran et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_chandran[["gex.rma"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_chandran <- addSlotMAE(mae_chandran, xcell = tmp)

# Friedrich et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_friedrich[["gex.logq"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_friedrich <- addSlotMAE(mae_friedrich, xcell = tmp)

# ICGCCA
tmp <- as.data.frame(immunedeconv::deconvolute(mae_icgcca[["gex.rma"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_icgcca <- addSlotMAE(mae_icgcca, xcell = tmp)

# Kunderfranco et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_kunderfranco[["gex.logr"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_kunderfranco <- addSlotMAE(mae_kunderfranco, xcell = tmp)

# Ren et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_ren[["gex.relz"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_ren <- addSlotMAE(mae_ren, xcell = tmp)

# Sun et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_sun[["gex.rma"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_sun <- addSlotMAE(mae_sun, xcell = tmp)

# Taylor et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_taylor[["gex.rma"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_taylor <- addSlotMAE(mae_taylor, xcell = tmp)

# TCGA
tmp <- as.data.frame(immunedeconv::deconvolute(mae_tcga[["gex.rsem.log"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_tcga <- addSlotMAE(mae_tcga, xcell = tmp)

# Wang et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_wang[["gex.rma"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
tmp <- as.matrix(tmp[,-1])
mae_wang <- addSlotMAE(mae_wang, xcell = tmp)

# Kim et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_kim[["gex.rma"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_kim <- addSlotMAE(mae_kim, xcell = tmp)

# Wallace et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_wallace[["gex.rma"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_wallace <- addSlotMAE(mae_wallace, xcell = tmp)

# Igc 
tmp <- as.data.frame(immunedeconv::deconvolute(mae_igc[["gex.rma"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_igc <- addSlotMAE(mae_igc, xcell = tmp)

# Weiner et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_weiner[["gex.rma"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_weiner <- addSlotMAE(mae_weiner, xcell = tmp)

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

# Barbieri et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_barbieri[["gex.relz"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_barbieri <- addSlotMAE(mae_barbieri, epic = tmp)

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

# Friedrich et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_friedrich[["gex.logq"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_friedrich <- addSlotMAE(mae_friedrich, epic = tmp)

# ICGCCA
tmp <- as.data.frame(immunedeconv::deconvolute(mae_icgcca[["gex.rma"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_icgcca <- addSlotMAE(mae_icgcca, epic = tmp)

# Kunderfranco et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_kunderfranco[["gex.logr"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_kunderfranco <- addSlotMAE(mae_kunderfranco, epic = tmp)

# Ren et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_ren[["gex.relz"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_ren <- addSlotMAE(mae_ren, epic = tmp)

# Sun et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_sun[["gex.rma"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_sun <- addSlotMAE(mae_sun, epic = tmp)

# Taylor et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_taylor[["gex.rma"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_taylor <- addSlotMAE(mae_taylor, epic = tmp)

# TCGA
tmp <- as.data.frame(immunedeconv::deconvolute(mae_tcga[["gex.rsem.log"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_tcga <- addSlotMAE(mae_tcga, epic = tmp)

# Wang et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_wang[["gex.rma"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_wang <- addSlotMAE(mae_wang, epic = tmp)

# Kim et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_kim[["gex.rma"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_kim <- addSlotMAE(mae_kim, epic = tmp)

# Wallace et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_wallace[["gex.rma"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_wallace <- addSlotMAE(mae_wallace, epic = tmp)

# IGC
tmp <- as.data.frame(immunedeconv::deconvolute(mae_igc[["gex.rma"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_igc <- addSlotMAE(mae_igc, epic = tmp)

# Wang et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_wang[["gex.rma"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_wang <- addSlotMAE(mae_wang, epic = tmp)

# Weiner et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_weiner[["gex.rma"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_weiner <- addSlotMAE(mae_weiner, epic = tmp)

# True et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_true[["gex.logr"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_true <- addSlotMAE(mae_true, epic = tmp)


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

# Barbieri et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_barbieri[["gex.relz"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_barbieri <- addSlotMAE(mae_barbieri, quantiseq = tmp)

# Barwick et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_barwick[["gex.logq"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_barwick <- addSlotMAE(mae_barwick, quantiseq = tmp)

# Chandran et al.  
tmp <- as.data.frame(immunedeconv::deconvolute(mae_chandran[["gex.rma"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_chandran <- addSlotMAE(mae_chandran, quantiseq = tmp)

# Friedrich et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_friedrich[["gex.logq"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_friedrich <- addSlotMAE(mae_friedrich, quantiseq = tmp)

# ICGCCA
tmp <- as.data.frame(immunedeconv::deconvolute(mae_icgcca[["gex.rma"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_icgcca <- addSlotMAE(mae_icgcca, quantiseq = tmp)

# Kunderfranco et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_kunderfranco[["gex.logr"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_kunderfranco <- addSlotMAE(mae_kunderfranco, quantiseq = tmp)

# Ren et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_ren[["gex.relz"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_ren <- addSlotMAE(mae_ren, quantiseq = tmp)

# Sun et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_sun[["gex.rma"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_sun <- addSlotMAE(mae_sun, quantiseq = tmp)

# Taylor et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_taylor[["gex.rma"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_taylor <- addSlotMAE(mae_taylor, quantiseq = tmp)

# TCGA
tmp <- as.data.frame(immunedeconv::deconvolute(mae_tcga[["gex.rsem.log"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_tcga <- addSlotMAE(mae_tcga, quantiseq = tmp)

# Wang et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_wang[["gex.rma"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_wang <- addSlotMAE(mae_wang, quantiseq = tmp)

# Kim et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_kim[["gex.rma"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_kim <- addSlotMAE(mae_kim, quantiseq = tmp)

# Wallace et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_wallace[["gex.rma"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_wallace <- addSlotMAE(mae_wallace, quantiseq = tmp)

# IGC
tmp <- as.data.frame(immunedeconv::deconvolute(mae_igc[["gex.rma"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_igc <- addSlotMAE(mae_igc, quantiseq = tmp)

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

# True et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_true[["gex.logr"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_true <- addSlotMAE(mae_true, quantiseq = tmp)



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

# Barbieri et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_barbieri[["gex.relz"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_barbieri <- addSlotMAE(mae_barbieri, mcp = tmp)

# Barwick et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_barwick[["gex.logq"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_barwick <- addSlotMAE(mae_barwick, mcp = tmp)

# Chandran et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_chandran[["gex.rma"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_chandran <- addSlotMAE(mae_chandran, mcp = tmp)

# Friedrich et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_friedrich[["gex.logq"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_friedrich <- addSlotMAE(mae_friedrich, mcp = tmp)

# ICGCCA
tmp <- as.data.frame(immunedeconv::deconvolute(mae_icgcca[["gex.rma"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_icgcca <- addSlotMAE(mae_icgcca, mcp = tmp)

# Kunderfranco et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_kunderfranco[["gex.logr"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_kunderfranco <- addSlotMAE(mae_kunderfranco, mcp = tmp)

# Ren et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_ren[["gex.relz"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_ren <- addSlotMAE(mae_ren, mcp = tmp)

# Sun et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_sun[["gex.rma"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_sun <- addSlotMAE(mae_sun, mcp = tmp)

# Taylor et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_taylor[["gex.rma"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_taylor <- addSlotMAE(mae_taylor, mcp = tmp)

# TCGA
tmp <- as.data.frame(immunedeconv::deconvolute(mae_tcga[["gex.rsem.log"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_tcga <- addSlotMAE(mae_tcga, mcp = tmp)

# Wang et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_wang[["gex.rma"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_wang <- addSlotMAE(mae_wang, mcp = tmp)

# Kim et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_kim[["gex.rma"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_kim <- addSlotMAE(mae_kim, mcp = tmp)

# Wallace et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_wallace[["gex.rma"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_wallace <- addSlotMAE(mae_wallace, mcp = tmp)

# IGC
tmp <- as.data.frame(immunedeconv::deconvolute(mae_igc[["gex.rma"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_igc <- addSlotMAE(mae_igc, mcp = tmp)

# Weiner et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_weiner[["gex.rma"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_weiner <- addSlotMAE(mae_weiner, mcp = tmp)

# True et al.
tmp <- as.data.frame(immunedeconv::deconvolute(mae_true[["gex.logr"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_true <- addSlotMAE(mae_true, mcp = tmp)



#####
##
## Genomic risk scores: Prolaris, OncotypeDX & Decipher
## AR scores as used by TCGA, originally presented in Hieronymus et al. 2006
##
#####


## Abida
mae_abida <- addSlotMAE(mae_abida,
	scores = rbind(
		decipher = curatedPCaData:::genomic_risk(mae_abida, slot = "gex.relz", test = "decipher", log=FALSE),
		oncotype = curatedPCaData:::genomic_risk(mae_abida, slot = "gex.relz", test = "oncotype", log=FALSE),
		prolaris = curatedPCaData:::genomic_risk(mae_abida, slot = "gex.relz", test = "prolaris", log=FALSE),
		ar_score = curatedPCaData:::genomic_score(mae_abida, slot = "gex.relz", test = "AR")
	)
)

## Barbieri
mae_barbieri <- addSlotMAE(mae_barbieri,
	scores = rbind(
		decipher = curatedPCaData:::genomic_risk(mae_barbieri, slot = "gex.relz", test = "decipher", log=FALSE),
		oncotype = curatedPCaData:::genomic_risk(mae_barbieri, slot = "gex.relz", test = "oncotype", log=FALSE),
		prolaris = curatedPCaData:::genomic_risk(mae_barbieri, slot = "gex.relz", test = "prolaris", log=FALSE),
		ar_score = curatedPCaData:::genomic_score(mae_barbieri, slot = "gex.relz", test = "AR")
	)
)

## Barwick
mae_barwick <- addSlotMAE(mae_barwick,
	scores = rbind(
		# Lack of gene overlap
		#decipher = curatedPCaData:::genomic_risk(mae_barwick, slot = "gex.logq", test = "decipher", log=FALSE),
		oncotype = curatedPCaData:::genomic_risk(mae_barwick, slot = "gex.logq", test = "oncotype", log=FALSE),
		prolaris = curatedPCaData:::genomic_risk(mae_barwick, slot = "gex.logq", test = "prolaris", log=FALSE),
		ar_score = curatedPCaData:::genomic_score(mae_barwick, slot = "gex.logq", test = "AR")
	)
)

## Chandran
mae_chandran <- addSlotMAE(mae_chandran,
	scores = rbind(
		decipher = curatedPCaData:::genomic_risk(mae_chandran, slot = "gex.rma", test = "decipher", log=FALSE),
		oncotype = curatedPCaData:::genomic_risk(mae_chandran, slot = "gex.rma", test = "oncotype", log=FALSE),
		prolaris = curatedPCaData:::genomic_risk(mae_chandran, slot = "gex.rma", test = "prolaris", log=FALSE),
		ar_score = curatedPCaData:::genomic_score(mae_chandran, slot = "gex.rma", test = "AR")
	)
)

## Friedrich
mae_friedrich <- addSlotMAE(mae_friedrich,
	scores = rbind(
		decipher = curatedPCaData:::genomic_risk(mae_friedrich, slot = "gex.logq", test = "decipher", log=FALSE),
		oncotype = curatedPCaData:::genomic_risk(mae_friedrich, slot = "gex.logq", test = "oncotype", log=FALSE),
		prolaris = curatedPCaData:::genomic_risk(mae_friedrich, slot = "gex.logq", test = "prolaris", log=FALSE),
		ar_score = curatedPCaData:::genomic_score(mae_friedrich, slot = "gex.logq", test = "AR")
	)
)

## ICGC-CA
mae_icgcca <- addSlotMAE(mae_icgcca,
	scores = rbind(
		decipher = curatedPCaData:::genomic_risk(mae_icgcca, slot = "gex.rma", test = "decipher", log=FALSE),
		oncotype = curatedPCaData:::genomic_risk(mae_icgcca, slot = "gex.rma", test = "oncotype", log=FALSE),
		prolaris = curatedPCaData:::genomic_risk(mae_icgcca, slot = "gex.rma", test = "prolaris", log=FALSE),
		ar_score = curatedPCaData:::genomic_score(mae_icgcca, slot = "gex.rma", test = "AR")
	)
)

## IGC
mae_igc <- addSlotMAE(mae_igc,
	scores = rbind(
		decipher = curatedPCaData:::genomic_risk(mae_igc, slot = "gex.rma", test = "decipher", log=FALSE),
		oncotype = curatedPCaData:::genomic_risk(mae_igc, slot = "gex.rma", test = "oncotype", log=FALSE),
		prolaris = curatedPCaData:::genomic_risk(mae_igc, slot = "gex.rma", test = "prolaris", log=FALSE),
		ar_score = curatedPCaData:::genomic_score(mae_igc, slot = "gex.rma", test = "AR")
	)
)

## Kim
mae_kim <- addSlotMAE(mae_kim,
	scores = rbind(
		decipher = curatedPCaData:::genomic_risk(mae_kim, slot = "gex.rma", test = "decipher", log=FALSE),
		oncotype = curatedPCaData:::genomic_risk(mae_kim, slot = "gex.rma", test = "oncotype", log=FALSE),
		prolaris = curatedPCaData:::genomic_risk(mae_kim, slot = "gex.rma", test = "prolaris", log=FALSE),
		ar_score = curatedPCaData:::genomic_score(mae_kim, slot = "gex.rma", test = "AR")
	)
)

## Kunderfranco
mae_kunderfranco <- addSlotMAE(mae_kunderfranco,
	scores = rbind(
		decipher = curatedPCaData:::genomic_risk(mae_kunderfranco, slot = "gex.logr", test = "decipher", log=FALSE),
		oncotype = curatedPCaData:::genomic_risk(mae_kunderfranco, slot = "gex.logr", test = "oncotype", log=FALSE),
		prolaris = curatedPCaData:::genomic_risk(mae_kunderfranco, slot = "gex.logr", test = "prolaris", log=FALSE),
		ar_score = curatedPCaData:::genomic_score(mae_kunderfranco, slot = "gex.logr", test = "AR")
	)
)

## Ren
mae_ren <- addSlotMAE(mae_ren,
	scores = rbind(
		decipher = curatedPCaData:::genomic_risk(mae_ren, slot = "gex.relz", test = "decipher", log=FALSE),
		oncotype = curatedPCaData:::genomic_risk(mae_ren, slot = "gex.relz", test = "oncotype", log=FALSE),
		prolaris = curatedPCaData:::genomic_risk(mae_ren, slot = "gex.relz", test = "prolaris", log=FALSE),
		ar_score = curatedPCaData:::genomic_score(mae_ren, slot = "gex.relz", test = "AR")
	)
)

## Sun
mae_sun <- addSlotMAE(mae_sun,
	scores = rbind(
		decipher = curatedPCaData:::genomic_risk(mae_sun, slot = "gex.rma", test = "decipher", log=FALSE),
		oncotype = curatedPCaData:::genomic_risk(mae_sun, slot = "gex.rma", test = "oncotype", log=FALSE),
		prolaris = curatedPCaData:::genomic_risk(mae_sun, slot = "gex.rma", test = "prolaris", log=FALSE),
		ar_score = curatedPCaData:::genomic_score(mae_sun, slot = "gex.rma", test = "AR")
	)
)

## Taylor
mae_taylor <- addSlotMAE(mae_taylor,
	scores = rbind(
		decipher = curatedPCaData:::genomic_risk(mae_taylor, slot = "gex.rma", test = "decipher", log=FALSE),
		oncotype = curatedPCaData:::genomic_risk(mae_taylor, slot = "gex.rma", test = "oncotype", log=FALSE),
		prolaris = curatedPCaData:::genomic_risk(mae_taylor, slot = "gex.rma", test = "prolaris", log=FALSE),
		ar_score = curatedPCaData:::genomic_score(mae_taylor, slot = "gex.rma", test = "AR")
	)
)

## TCGA 
mae_tcga <- addSlotMAE(mae_tcga,
	scores = rbind(
		decipher = curatedPCaData:::genomic_risk(mae_tcga, slot = "gex.rsem.log", test = "decipher", log=TRUE),
		oncotype = curatedPCaData:::genomic_risk(mae_tcga, slot = "gex.rsem.log", test = "oncotype", log=TRUE),
		prolaris = curatedPCaData:::genomic_risk(mae_tcga, slot = "gex.rsem.log", test = "prolaris", log=TRUE),
		ar_score = curatedPCaData:::genomic_score(mae_tcga, slot = "gex.rsem.log", test = "AR")
	)
)

## True
mae_true <- addSlotMAE(mae_true,
	scores = rbind(
		# Lack of gene overlap
		#decipher = curatedPCaData:::genomic_risk(mae_true, slot = "gex.logr", test = "decipher", log=FALSE),
		oncotype = curatedPCaData:::genomic_risk(mae_true, slot = "gex.logr", test = "oncotype", log=FALSE),
		prolaris = curatedPCaData:::genomic_risk(mae_true, slot = "gex.logr", test = "prolaris", log=FALSE),
		ar_score = curatedPCaData:::genomic_score(mae_true, slot = "gex.logr", test = "AR")
	)
)

## Wallace
mae_wallace <- addSlotMAE(mae_wallace,
	scores = rbind(
		decipher = curatedPCaData:::genomic_risk(mae_wallace, slot = "gex.rma", test = "decipher", log=FALSE),
		oncotype = curatedPCaData:::genomic_risk(mae_wallace, slot = "gex.rma", test = "oncotype", log=FALSE),
		prolaris = curatedPCaData:::genomic_risk(mae_wallace, slot = "gex.rma", test = "prolaris", log=FALSE),
		ar_score = curatedPCaData:::genomic_score(mae_wallace, slot = "gex.rma", test = "AR")
	)
)

## Wang
mae_wang <- addSlotMAE(mae_wang,
	scores = rbind(
		decipher = curatedPCaData:::genomic_risk(mae_wang, slot = "gex.rma", test = "decipher", log=FALSE),
		oncotype = curatedPCaData:::genomic_risk(mae_wang, slot = "gex.rma", test = "oncotype", log=FALSE),
		prolaris = curatedPCaData:::genomic_risk(mae_wang, slot = "gex.rma", test = "prolaris", log=FALSE),
		ar_score = curatedPCaData:::genomic_score(mae_wang, slot = "gex.rma", test = "AR")
	)
)

## Weiner
mae_weiner <- addSlotMAE(mae_weiner,
	scores = rbind(
		decipher = curatedPCaData:::genomic_risk(mae_weiner, slot = "gex.rma", test = "decipher", log=FALSE),
		oncotype = curatedPCaData:::genomic_risk(mae_weiner, slot = "gex.rma", test = "oncotype", log=FALSE),
		prolaris = curatedPCaData:::genomic_risk(mae_weiner, slot = "gex.rma", test = "prolaris", log=FALSE),
		ar_score = curatedPCaData:::genomic_score(mae_weiner, slot = "gex.rma", test = "AR")
	)
)



#####
##
## Purity estimate scores
## Potentially methods such as: DeMixT (tumor vs. ctrl), ISOpureR (tumor vs. ctrl), ESTIMATE (tumor-only)
##
#####

##
## For start, DeMixT for tumor-control design
##

# Remove DeMixT from current release
if(FALSE){
	## Chandran et al.
	chandran_demixt <- t(read.csv("data-raw/DeMixT_Chandran.txt", row.names=1))
	rownames(chandran_demixt) <- "demixt"
	mae_chandran <- addSlotMAE(mae_chandran, purity = chandran_demixt)

	## Friedrich et al.
	friedrich_demixt <- t(read.csv("data-raw/DeMixT_Friedrich.txt", row.names=1))
	rownames(friedrich_demixt) <- "demixt"
	mae_friedrich <- addSlotMAE(mae_friedrich, purity = friedrich_demixt)

	## Kunderfranco et al.
	kunderfranco_demixt <- t(read.csv("data-raw/DeMixT_Kunderfranco.txt", row.names=1))
	rownames(kunderfranco_demixt) <- "demixt"
	mae_kunderfranco <- addSlotMAE(mae_kunderfranco, purity = kunderfranco_demixt)

	## Taylor et al.
	taylor_demixt <- t(read.csv("data-raw/DeMixT_Taylor.txt", row.names=1))
	rownames(taylor_demixt) <- "demixt"
	mae_taylor <- addSlotMAE(mae_taylor, purity = taylor_demixt)

	## TCGA
	tcga_demixt <- t(read.csv("data-raw/DeMixT_TCGA.txt", row.names=1))
	rownames(tcga_demixt) <- "demixt"
	mae_tcga <- addSlotMAE(mae_tcga, purity = tcga_demixt)

	## Wallace et al.
	wallace_demixt <- t(read.csv("data-raw/DeMixT_Wallace.txt", row.names=1))
	rownames(wallace_demixt) <- "demixt"
	mae_wallace <- addSlotMAE(mae_wallace, purity = wallace_demixt)
}

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

# Re-save the processed datasets using the compression indicated at LazyDataCompression-field in DESCRIPTION
tools::resaveRdaFiles("data/mae_abida.rda", compress="xz")
tools::resaveRdaFiles("data/mae_barbieri.rda", compress="xz")
tools::resaveRdaFiles("data/mae_barwick.rda", compress="xz")
tools::resaveRdaFiles("data/mae_chandran.rda", compress="xz")
tools::resaveRdaFiles("data/mae_friedrich.rda", compress="xz")
tools::resaveRdaFiles("data/mae_icgcca.rda", compress="xz")
tools::resaveRdaFiles("data/mae_igc.rda", compress="xz")
tools::resaveRdaFiles("data/mae_kim.rda", compress="xz")
tools::resaveRdaFiles("data/mae_kunderfranco.rda", compress="xz")
tools::resaveRdaFiles("data/mae_ren.rda", compress="xz")
tools::resaveRdaFiles("data/mae_sun.rda", compress="xz")
tools::resaveRdaFiles("data/mae_taylor.rda", compress="xz")
tools::resaveRdaFiles("data/mae_tcga.rda", compress="xz")
tools::resaveRdaFiles("data/mae_true.rda", compress="xz")
tools::resaveRdaFiles("data/mae_wallace.rda", compress="xz")
tools::resaveRdaFiles("data/mae_wang.rda", compress="xz")
tools::resaveRdaFiles("data/mae_weiner.rda", compress="xz")

