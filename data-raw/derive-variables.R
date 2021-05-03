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


#####################################################
#####################################################
##                                                 ##
##                     XCELL                       ##
##                                                 ##
#####################################################
#####################################################

# Abida et al.

# 
tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_abida[["gex_capture"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_abida <- c(mae_abida, xcell_capture = tmp)

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_abida[["gex_polyA"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_abida <- c(mae_abida, xcell_polyA = tmp)

# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_abida, overwrite = TRUE)

# Barbieri et al.

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_barbieri[["gex"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_barbieri <- c(mae_barbieri, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_barbieri, overwrite = TRUE)

# Chandran et al.

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_chandran[["gex"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_chandran <- c(mae_chandran, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_chandran, overwrite = TRUE)

# Friedrich et al.

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_friedrich[["gex"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_friedrich <- c(mae_friedrich, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_friedrich, overwrite = TRUE)

# ICGCCA

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_icgcca[["gex"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_icgcca <- c(mae_icgcca, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_icgcca, overwrite = TRUE)

# Kunderfranco et al.

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_kunderfranco[["gex"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_kunderfranco <- c(mae_kunderfranco, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_kunderfranco, overwrite = TRUE)

# Ren et al.

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_ren[["gex"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_ren <- c(mae_ren, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_ren, overwrite = TRUE)

# Sun et al.

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_sun[["gex"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_sun <- c(curatedPCaData::mae_sun, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_sun, overwrite = TRUE)

# Taylor et al.

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_taylor[["gex"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_taylor <- c(curatedPCaData::mae_taylor, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_taylor, overwrite = TRUE)

# TCGA

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_tcga[["gex"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_tcga <- c(curatedPCaData::mae_tcga, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_tcga, overwrite = TRUE)

# Wallace et al.

#tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_wallace[["gex"]], method="xcell"))
#rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
#tmp <- as.matrix(tmp[,-1])
#mae_wallace <- c(mae_wallace, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_wallace, overwrite = TRUE)

# Weiner et al.

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_weiner[["gex"]], method="xcell"))
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

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_abida[["gex_capture"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_abida <- c(mae_abida, epic_capture = tmp)

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_abida[["gex_polyA"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_abida <- c(mae_abida, epic_polyA = tmp)

# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_abida, overwrite = TRUE)

# Barbieri et al.

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_barbieri[["gex"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_barbieri <- c(mae_barbieri, epic = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_barbieri, overwrite = TRUE)

# Chandran et al.

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_chandran[["gex"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_chandran <- c(mae_chandran, epic = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_chandran, overwrite = TRUE)

# Friedrich et al.

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_friedrich[["gex"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_friedrich <- c(mae_friedrich, epic = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_friedrich, overwrite = TRUE)

# ICGCCA

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_icgcca[["gex"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_icgcca <- c(mae_icgcca, epic = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_icgcca, overwrite = TRUE)

# Kunderfranco et al.

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_kunderfranco[["gex"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_kunderfranco <- c(mae_kunderfranco, epic = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_kunderfranco, overwrite = TRUE)

# Ren et al.

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_ren[["gex"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_ren <- c(mae_ren, epic = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_ren, overwrite = TRUE)

# Sun et al.

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_sun[["gex"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_sun <- c(curatedPCaData::mae_sun, epic = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_sun, overwrite = TRUE)

# Taylor et al.

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_taylor[["gex"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_taylor <- c(curatedPCaData::mae_taylor, epic = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_taylor, overwrite = TRUE)

# TCGA

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_tcga[["gex"]], method="epic"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_tcga <- c(curatedPCaData::mae_tcga, epic = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_tcga, overwrite = TRUE)

# Wallace et al.

#tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_wallace[["gex"]], method="epic"))
#rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
#tmp <- as.matrix(tmp[,-1])
#mae_wallace <- c(mae_wallace, epic = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_wallace, overwrite = TRUE)

# Weiner et al.

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_weiner[["gex"]], method="epic"))
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

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_abida[["gex_capture"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_abida <- c(mae_abida, quantiseq_capture = tmp)

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_abida[["gex_polyA"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_abida <- c(mae_abida, quantiseq_polyA = tmp)

# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_abida, overwrite = TRUE)

# Barbieri et al.

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_barbieri[["gex"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_barbieri <- c(mae_barbieri, quantiseq = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_barbieri, overwrite = TRUE)

# Chandran et al.  QUANTISEQ DOES NOT APPEAR TO WORK FOR CHANDRAN

#tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_chandran[["gex"]], method="quantiseq"))
#rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
#tmp <- as.matrix(tmp[,-1])
#mae_chandran <- c(mae_chandran, quantiseq = tmp)
# Save the derived new 'assay' types to the mae-object
#usethis::use_data(mae_chandran, overwrite = TRUE)

# Friedrich et al.

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_friedrich[["gex"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_friedrich <- c(mae_friedrich, quantiseq = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_friedrich, overwrite = TRUE)

# ICGCCA

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_icgcca[["gex"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_icgcca <- c(mae_icgcca, quantiseq = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_icgcca, overwrite = TRUE)

# Kunderfranco et al.

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_kunderfranco[["gex"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_kunderfranco <- c(mae_kunderfranco, quantiseq = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_kunderfranco, overwrite = TRUE)

# Ren et al.

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_ren[["gex"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_ren <- c(mae_ren, quantiseq = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_ren, overwrite = TRUE)

# Sun et al.

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_sun[["gex"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_sun <- c(curatedPCaData::mae_sun, quantiseq = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_sun, overwrite = TRUE)

# Taylor et al.

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_taylor[["gex"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_taylor <- c(curatedPCaData::mae_taylor, quantiseq = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_taylor, overwrite = TRUE)

# TCGA

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_tcga[["gex"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_tcga <- c(curatedPCaData::mae_tcga, quantiseq = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_tcga, overwrite = TRUE)

# Wallace et al.

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_wallace[["gex"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_wallace <- c(mae_wallace, quantiseq = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_wallace, overwrite = TRUE)

# Weiner et al.

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_weiner[["gex"]], method="quantiseq"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_weiner <- c(mae_weiner, quantiseq = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_weiner, overwrite = TRUE)



#####################################################
#####################################################
##                                                 ##
##                     MPC COUNTER                 ##
##                                                 ##
#####################################################
#####################################################

# Abida et al.

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_abida[["gex_capture"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_abida <- c(mae_abida, mcp_capture = tmp)

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_abida[["gex_polyA"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_abida <- c(mae_abida, mcp_polyA = tmp)

# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_abida, overwrite = TRUE)

# Barbieri et al.

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_barbieri[["gex"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_barbieri <- c(mae_barbieri, mcp = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_barbieri, overwrite = TRUE)

# Chandran et al.

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_chandran[["gex"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_chandran <- c(mae_chandran, mcp = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_chandran, overwrite = TRUE)

# Friedrich et al.

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_friedrich[["gex"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_friedrich <- c(mae_friedrich, mcp = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_friedrich, overwrite = TRUE)

# ICGCCA

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_icgcca[["gex"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_icgcca <- c(mae_icgcca, mcp = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_icgcca, overwrite = TRUE)

# Kunderfranco et al.

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_kunderfranco[["gex"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_wallace <- c(curatedPCaData::mae_wallace, mcp_counter = tmp)

mae_kunderfranco <- c(mae_kunderfranco, mcp = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_kunderfranco, overwrite = TRUE)

# Ren et al.

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_ren[["gex"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_ren <- c(mae_ren, mcp = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_ren, overwrite = TRUE)

# Sun et al.

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_sun[["gex"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_sun <- c(curatedPCaData::mae_sun, mcp = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_sun, overwrite = TRUE)

# Taylor et al.

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_taylor[["gex"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_taylor <- c(curatedPCaData::mae_taylor, mcp = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_taylor, overwrite = TRUE)

# TCGA

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_tcga[["gex"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_tcga <- c(curatedPCaData::mae_tcga, mcp = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_tcga, overwrite = TRUE)

# Wallace et al.

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_wallace[["gex"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_wallace <- c(mae_wallace, mcp = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_wallace, overwrite = TRUE)

# Weiner et al.

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_weiner[["gex"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
mae_weiner <- c(mae_weiner, mcp = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_weiner, overwrite = TRUE)


#####
##
## Genomic risk scores: Prolaris, OncotypeDX & Decipher
##
#####


rbind(
	Prolaris = curatedPCaData:::genomic_risk(curatedPCaData::mae_tcga, object = "gex", test = "Prolaris"),
	OncotypeDX = curatedPCaData:::genomic_risk(curatedPCaData::mae_tcga, object = "gex", test = "Oncotype DX"),
	Decipher = curatedPCaData:::genomic_risk(curatedPCaData::mae_tcga, object = "gex", test = "Decipher")
)

#####
##
## AR scores as used by TCGA, originally presented in Hieronymus et al. 2006
##
#####

