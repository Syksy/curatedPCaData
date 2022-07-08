#' Genomic risk scores
#'
#' Prolaris, Oncotype DX, Decipher transcriptome panels for prostate cancer risk
#'
#' @param mae MultiAssayExperiment object with gene expression available
#' @param slot Name of the Gene Expression slot, by default grepping for 'gex' prefix and picking the hit
#' @param test Type of test; available: "Prolaris", "Oncotype DX", and "Decipher" (case insensitive)
#' @param log Should data be log(x+1)-transformed prior to calculating the risk score
#'
#' @details https://bjui-journals.onlinelibrary.wiley.com/doi/10.1111/bju.14452 https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-14-690 https://www.nature.com/articles/s41391-019-0167-9
#'
#' @noRd
#' @keywords internal
genomic_risk <- function(mae, 
                         slot = grep("gex", names(mae), value=TRUE),
                         test = c("Prolaris", "Oncotype DX", "Decipher"),
                         log = FALSE # Should log-transformation be applied to the data for genomic risk score calculations; should not be utilized if data is transformed already
){
	# Internal function for automatically extracting the latest curatedPCaData::curatedPCaData_genes[,"Aliases"] for a specific hugo symbol
	expandAliases <- function(gene){
		if(length(which(curatedPCaData:::curatedPCaData_genes$hgnc_symbol == gene))>0){
			unique(c(gene, 
				unlist(
					strsplit(curatedPCaData:::curatedPCaData_genes[which(curatedPCaData:::curatedPCaData_genes$hgnc_symbol == gene),"Aliases"],";")[[1]]
				)
			))
		}else{
			gene
		}
	}
	# Internal function that extract expression value if it is present in matrix; otherwise imputes defined value, and gives a warning
	# List of aliases, and data matrix x with columns as genes
	extractGene <- function(genelist, x, impute = NA){
		if(any(genelist %in% colnames(x))){
			# Take the most recent alias hit
			x[,intersect(colnames(x), genelist)[1]]
		}else{
			warning(paste("Gene list", paste(genelist, collapse=";"), "not found from data matrix"))
			rep(impute, times=nrow(x))	
		}
	}

	# Prolaris
	prolaris_genes <- list(
		"FOXM1" = c("FOXM1"), 
		"CDC20" = c("CDC20"), 
		"CDKN3" = c("CDKN3"), 
		"CDC2" = c("CDC2"), 
		"KIF11" = c("KIF11"), 
		"KIAA0101" = c("KIAA0101"),
		"NUSAP1" = c("NUSAP1"), 
		"CENPF" = c("CENPF"), 
		"ASPM" = c("ASPM"), 
		"BUB1B" = c("BUB1B"), 
		"RRM2" = c("RRM2"), 
		"DLGAP5" = c("DLGAP5"),
		"BIRC5" = c("BIRC5"), 
		"KIF20A" = c("KIF20A"), 
		"PLK1" = c("PLK1"), 
		"TOP2A" = c("TOP2A"), 
		"TK1" = c("TK1"), 
		"PBK" = c("PBK"), 
		"ASF1B" = c("ASF1B"),
		"C18orf24" = c("C18orf24"), 
		"RAD54L" = c("RAD54L"), 
		"PTTG1" = c("PTTG1"), 
		"CDCA3" = c("CDCA3"), 
		"MCM10" = c("MCM10"), 
		"PRC1" = c("PRC1"),
		"DTL" = c("DTL"), 
		"CEP55" = c("CEP55"), 
		"RAD51" = c("RAD51"), 
		"CENPM" = c("CENPM"), 
		"CDCA8" = c("CDCA8"),
		"ORC6L" = c("ORC6L"), 
		"SKA1" = c("SKA1"),
		"ORC6" = c("ORC6"), 
		"CDK1" = c("CDK1")
	)
	prolaris_genes <- lapply(prolaris_genes, FUN=expandAliases)
	
  	# Oncotype DX
	oncotype_genes <- list(
		"AZGP1" = c("AZGP1"), 
		"KLK2" = c("KLK2"), 
		"SRD5A2" = c("SRD5A2"), 
		"FAM13C" = c("FAM13C"), 
		"FLNC" = c("FLNC"), 
		"GSN" = c("GSN"), 
		"TPM2" = c("TPM2"),
		"GSTM2" = c("GSTM2"), 
		"TPX2" = c("TPX2"), 
		"BGN" = c("BGN"), 
		"COL1A1" = c("COL1A1"), 
		"SFRP4" = c("SFRP4")
	)
	oncotype_genes <- lapply(oncotype_genes, FUN=expandAliases)

	# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3691249/
	# Table 2
	decipher_genes_over <- list(
		"CAMK2N1" = c("CAMK2N1"), 
		"EPPK1" = c("EPPK1"), 
		"IQGAP3" = c("IQGAP3"), 
		"LASP1" = c("LASP1"), 		
		"NFIB" = c("NFIB"), 
		"NUSAP1" = c("NUSAP1"),
		"PBX1" = c("PBX1"), 		
		"S1PR4" = c("S1PR4"), 		
		"THBS2" = c("THBS2"), 
		"UBE2C" = c("UBE2C"), 
		"ZWILCH" = c("ZWILCH")
	)
	decipher_genes_over <- lapply(decipher_genes_over, FUN=expandAliases)
	decipher_genes_under <- list(		
		"ANO7" = c("ANO7"), 
		"C6orf10" = c("C6orf10", "TSBP1"),
		"PCDH7" = c("PCDH7"), 
		"MYBPC1" = c("MYBPC1"), 
		"TSBP" = c("TSBP"), 
		"RABGAP1" = c("RABGAP1"), 
		"PCAT-32" = c("PCAT-32", "PCAT1"), 
		"PCAT-80" = c("PCAT-80", "GLYATL1P4"), 
		"TNFRSF19" = c("TNFRSF19")
	)
	decipher_genes_under <- lapply(decipher_genes_under, FUN=expandAliases)
	# Gene name annotations / changes from original publication
	# C6orf10 -> TSBP1
	# PCAT-32 -> PCAT1


	dat <- mae@ExperimentList[[slot]]
	dat <- t(dat) 
	dat <- as.data.frame(dat)
	# Log-transformation
	if(log){
		dat <- log2(dat+1)
	}
  	
  	# Prolaris
	if(base::tolower(test) == "prolaris"){
    
		if(length(intersect(colnames(dat), unlist(prolaris_genes))) != 31){
			warning("Prolaris risk score based off of ",
			      length(intersect(colnames(dat), unlist(prolaris_genes))),
			      " out of 31 genes")
		}
    
		risk <- dat[, colnames(dat) %in% unlist(prolaris_genes)]
		gene_med <- apply(risk, MARGIN=2, stats::median)
		#risk_centered <- risk - gene_med
		# Genes should be shifted in a row-wise manner
		risk_centered <- t(apply(risk, MARGIN=1, FUN=function(x) { x - gene_med }))
		
		# Squaring the median centered expression values 
		risk_centered <- risk_centered^2 

		risk_score <- apply(risk_centered, 1, mean)
		risk_score <- log2(risk_score)
		return(risk_score)
    
    	# Oncotype DX
	} else if (base::tolower(test) %in% c("oncotype dx", "oncotypedx", "oncotype")){
    
		if(!length(intersect(colnames(dat), unlist(oncotype_genes))) == 12){
			warning("The following required Oncotype DX genes: ", 
				paste(setdiff(names(oncotype_genes), colnames(dat)), collapse = ", "),
				" are not present in colnames of the data")
		}            
		#dat$TPX2_bounded <- ifelse(dat$TPX2 < 5, 5, dat$TPX2)
		#dat$SRD5A2_bounded <- ifelse(dat$SRD5A2 < 5.5, 5.5, dat$SRD5A2)

		# cellular_organization_module = dat$FLNC + dat$GSN + dat$TPM2 + dat$GSTM2
		# stromal_module = dat$BGN + dat$COL1A1 + dat$SFRP4
		# androgen_module = dat$FAM13C + dat$KLK2 + dat$SRD5A2_bounded + dat$AZGP1

		cellular_organization_module = 
			#(0.163*dat$FLNC) + 
			(0.163 * extractGene(genelist = oncotype_genes[["FLNC"]], x = dat, impute = 0)) + 
			#(0.504*dat$GSN) + 
			(0.504 * extractGene(genelist = oncotype_genes[["GSN"]], x = dat, impute = 0)) + 
			#(0.421*dat$TPM2) + 
			(0.421 * extractGene(genelist = oncotype_genes[["TPM2"]], x = dat, impute = 0)) + 
			#(0.394*dat$GSTM2)
			(0.394 * extractGene(genelist = oncotype_genes[["GSTM2"]], x = dat, impute = 0))
		stromal_module = 
			#(0.527*dat$BGN) + 
			(0.527 * extractGene(genelist = oncotype_genes[["BGN"]], x = dat, impute = 0)) + 
			#(0.457*dat$COL1A1) + 
			(0.457 * extractGene(genelist = oncotype_genes[["COL1A1"]], x = dat, impute = 0)) + 
			#(0.156*dat$SFRP4)
			(0.156 * extractGene(genelist = oncotype_genes[["SFRP4"]], x = dat, impute = 0))
		androgen_module = 
			#(0.634*dat$FAM13C) + 
			(0.634 * extractGene(genelist = oncotype_genes[["FAM13C"]], x = dat, impute = 0)) + 
			#(1.079*dat$KLK2) + 
			(1.079 * extractGene(genelist = oncotype_genes[["KLK2"]], x = dat, impute = 0)) + 
			#(0.997*dat$SRD5A2_bounded) + 
			# Lower bound 5.5
			(0.997 * unlist(lapply(extractGene(genelist = oncotype_genes[["SRD5A2"]], x = dat, impute = 0), FUN=function(x) {max(5.5, x) }))) + 
			#(0.642*dat$AZGP1)
			(0.642 * extractGene(genelist = oncotype_genes[["AZGP1"]], x = dat, impute = 0))
		proliferation_module = #dat$TPX2_bounded
			# Lower bound 5
			unlist(lapply(extractGene(genelist = oncotype_genes[["TPX2"]], x = dat, impute = 0), FUN=function(x) { max(5, x) }))

		risk_score = 0.735*stromal_module - 0.368*cellular_organization_module - 0.352*androgen_module + 0.095*proliferation_module
		names(risk_score) <- rownames(dat)	
			
		return(risk_score)
    
    	# Decipher
	} else if (base::tolower(test) %in% c("decipher", "decypher")) {
    
		if(length(intersect(colnames(dat), c(unlist(decipher_genes_over), unlist(decipher_genes_under)))) != 19){
			warning("The following required Decipher genes: ", 
				paste(setdiff(c(unlist(decipher_genes_over), unlist(decipher_genes_under)), colnames(dat)), collapse = ", "),
				" are not present in colnames of the GEX data")
		}
    
    		# Take intersection of available gene names and over and under expressed gene symbols and their aliases 
    		over <- intersect(colnames(dat), unlist(decipher_genes_over))
    		under <- intersect(colnames(dat), unlist(decipher_genes_under))
    
		# Intersect between Decipher risk score genes ideally available and current data matrix    
		risk <- dat[,intersect(colnames(dat),c(over, under))]

		# Median centering    
		gene_med <- apply(risk, MARGIN=2, stats::median)
		#risk_centered <- risk - gene_med
		# Genes should be shifted in a row-wise manner
		risk_centered <- t(apply(risk, MARGIN=1, FUN=function(x) { x - gene_med }))
    
		# average of the log2 normalized values for the 9 over-expressed targets
		c1 <- apply(risk_centered[,over], MARGIN=1, FUN=function(x) { mean(x, na.rm=TRUE) })
    
		# average of the log2 normalized values for the 9 under-expressed targets
		c2 <- apply(risk_centered[,under], MARGIN=1, FUN=function(x) { mean(x, na.rm=TRUE) })

		# Risk score as the difference between over- and underrepresented genes    
		risk_score <- c1-c2
    
		return(risk_score)
    
	} else {
		stop(paste("Invalid genomic risk score name:", test))
	}
}


#' Various genomic scores
#'
#' AR score by Hieronymus et al 2006 as used by TCGA 2015
#' Quote: "To address these questions, we sought to infer the AR output of tumors 
#' by calculating an AR activity score from the expression pattern of 20 genes that 
#' are experimentally validated AR transcriptional targets (Hieronymus et al., 2006)."
#'
#' @noRd
#' @keywords internal
genomic_score <- function(
			mae, # MultiAssayExperiment object
			slot = "gex", # Slot inside MAE object to use as the gene expression
			test = "AR", # Test/score to calculate; by default Androgen Receptor (AR) score is calculated
			verbose = TRUE # Warnings for not found symbols etc
){
	# TCGA methodology for AR output score analysis: (Section 6 in supplementary of https://www.cell.com/cms/10.1016/j.cell.2015.10.025/attachment/70a60372-cdaf-4c72-aa6d-ded4b33ce5a0/mmc1.pdf )
	# "The AR output score is derived from the mRNA expression of genes that are experimentally
	# validated AR transcriptional targets (Hieronymus et al., 2006). Precisely, a list of 20 genes
	# upregulated in LNCaP cells stimulated with the synthetic androgen R1881 was used as a gene
	# signature of androgen-induced genes. An AR output score was defined by the quantification of
	# the composite expression of this 20-gene signature in each sample. Here, we measured
	# differential AR activity between genomic subtypes (ERG, ETV1/4/FLI1, SPOP, FOXA1, other,
	# normal prostate). To this aim, we computed a Z-score for the expression of each gene in each
	# sample by subtracting the pooled mean from the RNA-seq expression values and dividing by
	# the pooled standard deviation."
	

	# https://www.sciencedirect.com/science/article/pii/S1535610806002820
	# Fig 1B
	# " A gene expression signature of androgen stimulation was defined from gene expression profiles of LNCaP cells stimulated with the synthetic androgen R1881 for 12 hr and 24 hr, 
	# as compared to androgen-deprived LNCaP cells. The 27 gene signature contains both androgen-induced and androgen-repressed genes, shown here by row-normalized heat map."
	#
	## Genes as they were in Hieronymus et al 2006
	if(FALSE){
		hieronymus_genes_up <- c(
			"KLK3", #"PSA", # PSA -> KLK3 gene
			"TMPRSS2", "NKX3.1", # "NKX3-1", # Aliases
			"KLK2", "GNMT", "PMEPA1", # "TMEPAI", # Updated annotation; TMEPAI -> PMEPA1
			"MPHOSHP9", #"MPHOS9", # MPHOS9 -> MPHOSPH9
			"ZBTB10", "EAF2", 
			"CENPN", # "BM039", # BM039 -> CENPN
			#"SARG", # SARG not found; could be C1orf116
			"ACSL3", "PTGER4", "ABCC4",
			"NNMT", "ADAM7", "FKBP5", "ELL2", "MED28", "HERC3", "MAF")
		# Based on Fib 1C ELL2 might fit better into down than up
		hieronymus_genes_dn <- c("TNK1", "GLRA2", "MAPRE2", "PIP5K2B", "MAN1A1", "CD200")
	}
	# TCGA version of Hieronymus AR-genes (panel A rows): https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4695400/figure/F4/ 
	if(FALSE){
		hieronymus_genes_up <- c(
			"KLK3", "KLK2", "PMEPA1",
			"ABCC4", 
			"NKX3-1", "NKX3.1", # Two naming conventions, '-' replaced with '.'
			"C1orf116", # Not conventionally found from GEX
			"FKBP5", "ACSL3", "ZBTB10", "HERC3", 
			"PTGER4", "MPHOSPH9", "EAF2", "MED28", "NNMT", "MAF",
			"GNMT", "CENPN", "ELL2", "TMPRSS2"
		)	
	}
	# TCGA version of AR score supporting gene aliases and different naming conventions
	
	if(base::tolower(test) %in%  c("ar", "ar score", "ar-score", "ar_score")){		
		# Aliases queried using https://www.genecards.org/
		ar_genes <- list(
			"KLK3" = c("KLK3", "PSA", "APS", "KLK2A1"), # Possibly HK3; ambiguous
			"KLK2" = c("KLK2", "HGK-1", "HGK.1", "KLK2A2"), # Possibly HK2; ambiguous
			"PMEPA1" = c("PMEPA1", "STAG1", "TMEPAI"),
			"ABCC4" = c("ABCC4", "MRP4", "MOATB", "MOAT-B", "MOAT.B"),
			"NKX3-1" = c("NKX3-1", "NKX3.1", "BAPX2", "NKX3A"),
			"C1orf116" = c("C1orf116", "SARG", "FLJ36507", "MGC2742", "MGC4309"),
			"FKBP5" = c("FKBP5", "FKBP51", "FKBP54", "FKBP-51", "FKBP.51", "AIG6", "FKBP-5", "FKBP.5"),
			"ACSL3" = c("ACSL3", "FACL3", "LACS3"),
			"ZBTB10" = c("ZBTB10", "RINZFC", "RINZF"),
			"HERC3" = c("HERC3"),
			"PTGER4" = c("PTGER4", "EP4", "P4R"), # Possibly PTGER2; ambiguous
			"MPHOSPH9" = c("MPHOSPH9", "MPHOS9", "MPP9", "MPP-9", "MPP.9"),
			"EAF2" = c("EAF2", "TRAITS", "BM040", "U19"),
			"MED28" = c("MED28", "EG1"),
			"NNMT" = c("NNMT"),
			"MAF" = c("MAF", "CTRCT21", "AYGRP", "CCA4", "C-MAF", "C.MAF"),
			"GNMT" = c("GNMT"),
			"CENPN" = c("CENPN", "ICEN32", "BM039", "FLJ13607", "FLJ22660"),
			"ELL2" = c("ELL2", "MRCCAT1"),
			"TMPRSS2" = c("TMPRSS2", "PRSS10")
		)
	
		missing <- 0

		# Pooled normalization - pooling within sample over all genes as done in TCGA
		gex <- mae[[slot]]
		gez <- t(apply(gex, MARGIN=1, FUN=function(z) { 
			scale(z, center=TRUE, scale=TRUE)
		}))
		dimnames(gez) <- dimnames(gex)
		
		# TCGA computed AR scores based on just up regulated genes; sum of z-scores from pooled normalization
		# Check for gene name overlaps and warn of missing ones
		if(verbose){
			lapply(1:length(ar_genes), FUN=function(z){
				if(length(intersect(rownames(gez), ar_genes[[z]]))>1){
					warning(paste("Warning! More than one gene name alias found for:", names(ar_genes)[z]))
				}
				if(length(intersect(rownames(gez), ar_genes[[z]]))==0){
					warning(paste("Warning! No gene name from alias list found for:", names(ar_genes)[z]))
				}
			})
		}
		# lapply over patients
		res <- unlist(lapply(colnames(gez), FUN=function(patient){
			# lapply over genes for AR score and their aliases
			sum(unlist(lapply(ar_genes, FUN=function(z){
				gez[intersect(rownames(gez), z), patient]
			})), na.rm=TRUE)
		}))
		names(res) <- colnames(gex)
		res		
	}else{
		stop(paste("Unknown genomic score parameter:", type))
	}
}

