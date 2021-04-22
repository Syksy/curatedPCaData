#####
##
## -- DERIVING NEW VARIABLES --
## Set of scripts that derive new variables from the original raw data
##
#####

library(immunedeconv)

#####
##
## xCell - immmune cell decomposition (rank-base approach; rather universal)
##
#####

##
# TCGA
##

# Example, xCell based on gene expression, and then concatenate this to mae_tcga

# Transform from tibble to more universal formats such as a data.frame or a matrix
tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_tcga[["gex"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_tcga <- c(curatedPCaData::mae_tcga, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_tcga, overwrite = TRUE)

##
# Taylor et al.
##

# Transform from tibble to more universal formats such as a data.frame or a matrix
tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_taylor[["gex"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_taylor <- c(curatedPCaData::mae_taylor, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_taylor, overwrite = TRUE)

##
# Sun et al.
##

# Transform from tibble to more universal formats such as a data.frame or a matrix
tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_sun[["gex"]], method="xcell"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_sun <- c(curatedPCaData::mae_sun, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_sun, overwrite = TRUE)

##
# Friedrich et al.
##

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_friedrich[["gex"]], method="xcell"))

rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])

mae_friedrich <- c(mae_friedrich, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_friedrich, overwrite = TRUE)

##
# Wallace et al.
##

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_wallace[["gex"]], method="xcell"))

rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])

mae_wallace <- c(mae_wallace, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
<<<<<<< HEAD
## !!
## !! TDL: NOTE! There was a typo below (wallce instead of wallace) !!
## !!
usethis::use_data(mae_wallace, overwrite = TRUE)
=======
usethis::use_data(mae_wallace, overwrite = TRUE)

###
## Weiner et al.
###

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_weiner[["gex"]], method="xcell"))

rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])

mae_weiner <- c(mae_weiner, xcell = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_weiner, overwrite = TRUE)
>>>>>>> origin/f_branch



##
# -- Recompile and re-install package, so expanded mae-objects are up to date --
##



#####
##
## MCP-counter (sums of gene expression values; not sensitive to normalization method)
##
####

##
# TCGA
##

# Example, xCell based on gene expression, and then concatenate this to mae_tcga

# Transform from tibble to more universal formats such as a data.frame or a matrix
tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_tcga[["gex"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_tcga <- c(curatedPCaData::mae_tcga, mcp_counter = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_tcga, overwrite = TRUE)

##
# Taylor et al.
##

# Transform from tibble to more universal formats such as a data.frame or a matrix
tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_taylor[["gex"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_taylor <- c(curatedPCaData::mae_taylor, mcp_counter = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_taylor, overwrite = TRUE)

##
# Sun et al.
##

# Transform from tibble to more universal formats such as a data.frame or a matrix
tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_sun[["gex"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_sun <- c(curatedPCaData::mae_sun, mcp_counter = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_sun, overwrite = TRUE)

##
# Friedrich et al.
##

# Transform from tibble to more universal formats such as a data.frame or a matrix
tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_friedrich[["gex"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
mae_friedrich <- c(curatedPCaData::mae_friedrich, mcp_counter = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_friedrich, overwrite = TRUE)




##
# Wallace et al.
##

# Transform from tibble to more universal formats such as a data.frame or a matrix
tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_wallace[["gex"]], method="mcp_counter"))
rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])
# Concatenate the new results to the MAE-object
# !!
# !! TDL: NOTE! There was a typo below (wallce instead of wallace) !!
# !!
mae_wallace <- c(curatedPCaData::mae_wallace, mcp_counter = tmp)
# Save the derived new 'assay' types to the mae-object
usethis::use_data(mae_wallace, overwrite = TRUE)


##
# Weiner et al.
##

tmp <- as.data.frame(immunedeconv::deconvolute(curatedPCaData::mae_weiner[["gex"]], method="mcp_counter"))

rownames(tmp) <- tmp$cell_type
# Omit cell type column and store only data of cell type populations
tmp <- as.matrix(tmp[,-1])

mae_weiner <- c(mae_weiner, mcp_counter = tmp)
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

