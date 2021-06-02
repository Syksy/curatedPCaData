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
cibersort_kunderfranco <- t(cibersort_kunderfranco)
colnames(cibersort_kunderfranco)<- cibersort_kunderfranco[1,]
cibersort_kunderfranco<-cibersort_kunderfranco[-1,]
cibersort_kunderfranco<-as.matrix(cibersort_kunderfranco)
save(cibersort_kunderfranco, file="data-raw/cibersort_kunderfranco.RData")
#mae_kunderfranco <- c(mae_kunderfranco, cibersort = cibersort_kunderfranco)
mae_kunderfranco <- create_mae(study_name = "kunderfranco")
usethis::use_data(mae_kunderfranco, overwrite = TRUE,internal=TRUE)

# Sun et al.

cibersort_sun<-rio::import("data-raw/CIBERSORTx_sun_Results.csv")
cibersort_sun <- t(cibersort_sun)
colnames(cibersort_sun)<- cibersort_sun[1,]
cibersort_sun<-cibersort_sun[-1,]
cibersort_sun<-as.matrix(cibersort_sun)
save(cibersort_sun, file="data-raw/cibersort_sun.RData")
#mae_sun <- c(mae_sun, cibersort = cibersort_sun)
mae_sun <- create_mae(study_name = "sun")
usethis::use_data(mae_sun, overwrite = TRUE,internal=TRUE)

# Wang et al.

cibersort_wang<-rio::import("data-raw/CIBERSORTx_wang_Results.csv")
cibersort_wang <- t(cibersort_wang)
colnames(cibersort_wang)<- cibersort_wang[1,]
cibersort_wang<-cibersort_wang[-1,]
cibersort_wang<-as.matrix(cibersort_wang)
save(cibersort_wang, file="data-raw/cibersort_wang.RData")
#mae_wang <- c(mae_wang, cibersort = cibersort_wang)
mae_wang <- create_mae(study_name = "wang")
usethis::use_data(mae_wang, overwrite = TRUE,internal=TRUE)

# Kim et al.

cibersort_kim<-rio::import("data-raw/CIBERSORTx_kim_Results.csv")
cibersort_kim <- t(cibersort_kim)
colnames(cibersort_kim)<- cibersort_kim[1,]
cibersort_kim<-cibersort_kim[-1,]
cibersort_kim<-as.matrix(cibersort_kim)
save(cibersort_kim, file="data-raw/cibersort_kim.RData")
#mae_kim <- c(mae_kim, cibersort = cibersort_kim)
mae_kim <- create_mae(study_name = "kim")
usethis::use_data(mae_kim, overwrite = TRUE,internal=TRUE)

# Barbieri et al.

cibersort_barbieri<-rio::import("data-raw/CIBERSORTx_barbieri_Results.csv")
cibersort_barbieri <- t(cibersort_barbieri)
colnames(cibersort_barbieri)<- cibersort_barbieri[1,]
cibersort_barbieri<-cibersort_barbieri[-1,]
cibersort_barbieri<-as.matrix(cibersort_barbieri)
save(cibersort_barbieri, file="data-raw/cibersort_barbieri.RData")
#mae_barbieri <- c(mae_barbieri, cibersort = cibersort_barbieri)
mae_barbieri <- create_mae(study_name = "barbieri")
usethis::use_data(mae_barbieri, overwrite = TRUE,internal=TRUE)



library(immunedeconv) # NECESSARY!

#####################################################
#####################################################
##                                                 ##
##                     XCELL                       ##
##                                                 ##
#####################################################
#####################################################


# Abida et al.

# 
tmp <- as.data.frame(immunedeconv::deconvolute(mae_abida[["gex_capture"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_abida <- c(mae_abida, xcell_capture = tmp)

tmp <- as.data.frame(immunedeconv::deconvolute(mae_abida[["gex_polyA"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_abida <- c(mae_abida, xcell_polyA = tmp)

# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_abida, overwrite = TRUE)

# Barbieri et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_barbieri[["gex"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_barbieri <- c(mae_barbieri, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_barbieri, overwrite = TRUE,internal = TRUE)

# Chandran et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_chandran[["gex"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_chandran <- c(mae_chandran, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_chandran, overwrite = TRUE)

# Friedrich et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_friedrich[["gex"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_friedrich <- c(mae_friedrich, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_friedrich, overwrite = TRUE)

# ICGCCA

tmp <- as.data.frame(immunedeconv::deconvolute(mae_icgcca[["gex"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_icgcca <- c(mae_icgcca, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_icgcca, overwrite = TRUE)

# Kunderfranco et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_kunderfranco[["gex"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_kunderfranco <- c(mae_kunderfranco, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_kunderfranco, overwrite = TRUE,internal = TRUE)

# Ren et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_ren[["gex"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_ren <- c(mae_ren, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_ren, overwrite = TRUE)

# Sun et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_sun[["gex"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_sun <- c(curatedPCaData::mae_sun, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_sun, overwrite = TRUE)

# Taylor et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_taylor[["gex"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_taylor <- c(curatedPCaData::mae_taylor, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_taylor, overwrite = TRUE,internal = TRUE)

# TCGA

tmp <- as.data.frame(immunedeconv::deconvolute(mae_tcga[["gex"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_tcga <- c(curatedPCaData::mae_tcga, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_tcga, overwrite = TRUE,internal = TRUE)

# Wang et al.

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_wang[["gex"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_wang <- c(curatedPCaData::mae_wang, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_wang, overwrite = TRUE,internal = TRUE)

# Kim et al.

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_kim[["gex"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_kim <- c(curatedPCaData::mae_kim, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_kim, overwrite = TRUE,internal = TRUE)

# Wallace et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_wallace[["gex"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_wallace <- c(mae_wallace, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_wallace, overwrite = TRUE)

# Weiner et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_weiner[["gex"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_weiner <- c(mae_weiner, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_weiner, overwrite = TRUE)

#####################################################
#####################################################
##                                                 ##
##                     EPIC                        ##
##                                                 ##
#####################################################
#####################################################

# Abida et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_abida[["gex_capture"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_abida <- c(mae_abida, epic_capture = tmp)

tmp <- as.data.frame(immunedeconv::deconvolute(mae_abida[["gex_polyA"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_abida <- c(mae_abida, epic_polyA = tmp)

# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_abida, overwrite = TRUE)

# Barbieri et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_barbieri[["gex"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_barbieri <- c(mae_barbieri, epic = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_barbieri, overwrite = TRUE,internal = TRUE)

# Chandran et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_chandran[["gex"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_chandran <- c(mae_chandran, epic = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_chandran, overwrite = TRUE)

# Friedrich et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_friedrich[["gex"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_friedrich <- c(mae_friedrich, epic = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_friedrich, overwrite = TRUE)

# ICGCCA

tmp <- as.data.frame(immunedeconv::deconvolute(mae_icgcca[["gex"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_icgcca <- c(mae_icgcca, epic = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_icgcca, overwrite = TRUE)

# Kunderfranco et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_kunderfranco[["gex"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_kunderfranco <- c(mae_kunderfranco, epic = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_kunderfranco, overwrite = TRUE,internal = TRUE)


# Ren et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_ren[["gex"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_ren <- c(mae_ren, epic = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_ren, overwrite = TRUE)

# Sun et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_sun[["gex"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_sun <- c(curatedPCaData::mae_sun, epic = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_sun, overwrite = TRUE)

# Taylor et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_taylor[["gex"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_taylor <- c(curatedPCaData::mae_taylor, epic = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_taylor, overwrite = TRUE,internal = TRUE)

# TCGA

tmp <- as.data.frame(immunedeconv::deconvolute(mae_tcga[["gex"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_tcga <- c(curatedPCaData::mae_tcga, epic = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_tcga, overwrite = TRUE,internal = TRUE)

# Wang et al.

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_wang[["gex"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_wang <- c(curatedPCaData::mae_wang, epic = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_wang, overwrite = TRUE,internal = TRUE)

# Kim et al.

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_kim[["gex"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_kim <- c(curatedPCaData::mae_kim, epic = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_kim, overwrite = TRUE,internal = TRUE)

# Wallace et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_wallace[["gex"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_wallace <- c(mae_wallace, epic = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_wallace, overwrite = TRUE)

# Weiner et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_weiner[["gex"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_weiner <- c(mae_weiner, epic = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_weiner, overwrite = TRUE)


#####################################################
#####################################################
##                                                 ##
##                     QUANTISEQ                   ##
##                                                 ##
#####################################################
#####################################################

# Abida et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_abida[["gex_capture"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_abida <- c(mae_abida, quantiseq_capture = tmp)

tmp <- as.data.frame(immunedeconv::deconvolute(mae_abida[["gex_polyA"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_abida <- c(mae_abida, quantiseq_polyA = tmp)

# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_abida, overwrite = TRUE)

# Barbieri et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_barbieri[["gex"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_barbieri <- c(mae_barbieri, quantiseq = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_barbieri, overwrite = TRUE,internal = TRUE)

# Chandran et al.  QUANTISEQ DOES NOT APPEAR TO WORK FOR CHANDRAN

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_chandran[["gex"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_chandran <- c(mae_chandran, quantiseq = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_chandran, overwrite = TRUE)

# Friedrich et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_friedrich[["gex"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_friedrich <- c(mae_friedrich, quantiseq = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_friedrich, overwrite = TRUE)

# ICGCCA

tmp <- as.data.frame(immunedeconv::deconvolute(mae_icgcca[["gex"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_icgcca <- c(mae_icgcca, quantiseq = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_icgcca, overwrite = TRUE)

# Kunderfranco et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_kunderfranco[["gex"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_kunderfranco <- c(mae_kunderfranco, quantiseq = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_kunderfranco, overwrite = TRUE,internal=TRUE)

# Ren et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_ren[["gex"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_ren <- c(mae_ren, quantiseq = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_ren, overwrite = TRUE)

# Sun et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_sun[["gex"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_sun <- c(curatedPCaData::mae_sun, quantiseq = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_sun, overwrite = TRUE,internal = TRUE)

# Taylor et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_taylor[["gex"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_taylor <- c(curatedPCaData::mae_taylor, quantiseq = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_taylor, overwrite = TRUE,internal = TRUE)

# TCGA

tmp <- as.data.frame(immunedeconv::deconvolute(mae_tcga[["gex"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_tcga <- c(curatedPCaData::mae_tcga, quantiseq = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_tcga, overwrite = TRUE,internal=TRUE)

# Wang et al.

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_wang[["gex"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_wang <- c(curatedPCaData::mae_wang, quantiseq = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_wang, overwrite = TRUE,internal=TRUE)

# Kim et al.

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_kim[["gex"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_kim <- c(curatedPCaData::mae_kim, quantiseq = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_kim, overwrite = TRUE,internal=TRUE)

# Wallace et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_wallace[["gex"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_wallace <- c(mae_wallace, quantiseq = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_wallace, overwrite = TRUE)

# Weiner et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_weiner[["gex"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_weiner <- c(mae_weiner, quantiseq = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_weiner, overwrite = TRUE)



#####################################################
#####################################################
##                                                 ##
##                     MCP COUNTER                 ##
##                                                 ##
#####################################################
#####################################################

# Abida et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_abida[["gex_capture"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_abida <- c(mae_abida, mcp_capture = tmp)

tmp <- as.data.frame(immunedeconv::deconvolute(mae_abida[["gex_polyA"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_abida <- c(mae_abida, mcp_polyA = tmp)

# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_abida, overwrite = TRUE)

# Barbieri et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_barbieri[["gex"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_barbieri <- c(mae_barbieri, mcp = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_barbieri, overwrite = TRUE,internal = TRUE)

# Chandran et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_chandran[["gex"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_chandran <- c(mae_chandran, mcp = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_chandran, overwrite = TRUE)

# Friedrich et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_friedrich[["gex"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_friedrich <- c(mae_friedrich, mcp = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_friedrich, overwrite = TRUE)

# ICGCCA

tmp <- as.data.frame(immunedeconv::deconvolute(mae_icgcca[["gex"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_icgcca <- c(mae_icgcca, mcp = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_icgcca, overwrite = TRUE)

# Kunderfranco et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_kunderfranco[["gex"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
#mae_wallace <- c(curatedPCaData::mae_wallace, mcp_counter = tmp)
mae_kunderfranco <- c(curatedPCaData::mae_kunderfranco, mcp = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_kunderfranco, overwrite = TRUE,internal=TRUE)

# Ren et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_ren[["gex"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_ren <- c(mae_ren, mcp = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_ren, overwrite = TRUE)

# Sun et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_sun[["gex"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_sun <- c(curatedPCaData::mae_sun, mcp = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_sun, overwrite = TRUE,internal=TRUE)

# Taylor et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_taylor[["gex"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_taylor <- c(curatedPCaData::mae_taylor, mcp = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_taylor, overwrite = TRUE,internal = TRUE)

# TCGA

tmp <- as.data.frame(immunedeconv::deconvolute(mae_tcga[["gex"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_tcga <- c(curatedPCaData::mae_tcga, mcp = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_tcga, overwrite = TRUE,internal = TRUE)

# Wang et al.

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_wang[["gex"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_wang <- c(curatedPCaData::mae_wang, mcp = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_wang, overwrite = TRUE,internal = TRUE)

# Kim et al.

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_kim[["gex"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_kim <- c(curatedPCaData::mae_kim, mcp = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_kim, overwrite = TRUE,internal = TRUE)

# Wallace et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_wallace[["gex"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_wallace <- c(mae_wallace, mcp = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_wallace, overwrite = TRUE)

# Weiner et al.

tmp <- as.data.frame(immunedeconv::deconvolute(mae_weiner[["gex"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_weiner <- c(mae_weiner, mcp = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_weiner, overwrite = TRUE)


#####
##
## Genomic risk scores: Prolaris, OncotypeDX & Decipher
## AR scores as used by TCGA, originally presented in Hieronymus et al. 2006
#####


## Abida

mae_abida <- c(mae_abida,
	scores = rbind(
		Prolaris = curatedPCaData:::genomic_risk(curatedPCaData::mae_abida, object = "gex_capture", test = "Prolaris"),
		AR_score = curatedPCaData:::genomic_score(curatedPCaData::mae_abida, object = "gex_capture", test = "AR")
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
usethis::use_data(mae_tcga, overwrite = TRUE)

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
