#' Genomic risk scores
#'
#' Prolaris, Oncotype DX, Decipher genomic panels for prostate cancer risk
#'
#' @details https://bjui-journals.onlinelibrary.wiley.com/doi/10.1111/bju.14452 https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-14-690 https://www.nature.com/articles/s41391-019-0167-9
#'
genomic_risk <- function(mae, 
                         slot = "gex",
                         test = c("Prolaris", "Oncotype DX", "Decipher"),
                         log = TRUE # Should log-transformation be applied to the data for genomic risk score calculations; should not be utilized if data transformed already
){
  
  prolaris_genes <- c("FOXM1", "CDC20", "CDKN3", "CDC2", "KIF11", "KIAA0101",
                      "NUSAP1", "CENPF", "ASPM", "BUB1B", "RRM2", "DLGAP5",
                      "BIRC5", "KIF20A", "PLK1", "TOP2A", "TK1", "PBK", "ASF1B",
                      "C18orf24", "RAD54L", "PTTG1", "CDCA3", "MCM10", "PRC1",
                      "DTL", "CEP55", "RAD51", "CENPM", "CDCA8", "ORC6L", "SKA1",
                      "ORC6", "CDK1")
  
  oncotype_genes <- c("AZGP1", "KLK2", "SRD5A2", "FAM13C", "FLNC", "GSN", "TPM2",
                      "GSTM2", "TPX2", "BGN", "COL1A1", "SFRP4")
  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3691249/
  # Table 2
  decipher_genes <- c("LASP1", "IQGAP3", "NFIB", "S1PR4", "THBS2", "ANO7", 
                      "PCDH7", "MYBPC1", "EPPK1", "TSBP", "PBX1", "NUSAP1",
                      "ZWILCH", "UBE2C", "CAMK2N1", "RABGAP1", "PCAT-32", 
                      "PCAT-80", "TNFRSF19", "C6orf10")
  # Gene name annotations / changes from original publication
  # C6orf10 -> TSBP1
  # PCAT-32 -> PCAT1
  
  
  dat <- mae@ExperimentList[[slot]]
  dat <- t(dat) 
  dat <- as.data.frame(dat)
  
  if(base::tolower(test) == "prolaris"){
    
    if(length(intersect(colnames(dat), prolaris_genes)) != 31){
      warning("Prolaris risk score based off of ",
              length(intersect(colnames(dat), prolaris_genes)),
              " out of 31 genes")
    }
    
    risk <- dat[, colnames(dat) %in% prolaris_genes]
    gene_med <- apply(risk, 2, stats::median)
    risk_centered <- risk - gene_med
    risk_centered <- risk_centered^2 # squaring the median centered expression values 
    
    risk_score <- apply(risk_centered, 1, mean)
    risk_score <- log2(risk_score)
    return(risk_score)
    
  } else if (base::tolower(test) %in% c("oncotype dx", "oncotypedx", "oncotype")){
    
    if(length(intersect(colnames(dat), oncotype_genes)) == 12){
      
      if(log){
      	dat <- log2(dat)
      }
      
      dat$TPX2_bounded <- ifelse(dat$TPX2 < 5, 5, dat$TPX2)
      dat$SRD5A2_bounded <- ifelse(dat$SRD5A2 < 5.5, 5.5, dat$SRD5A2)
      
      # cellular_organization_module = dat$FLNC + dat$GSN + dat$TPM2 + dat$GSTM2
      # stromal_module = dat$BGN + dat$COL1A1 + dat$SFRP4
      # androgen_module = dat$FAM13C + dat$KLK2 + dat$SRD5A2_bounded + dat$AZGP1
      
      cellular_organization_module = (0.163*dat$FLNC) + 
        (0.504*dat$GSN) + (0.421*dat$TPM2) + (0.394*dat$GSTM2)
      stromal_module = (0.527*dat$BGN) + (0.457*dat$COL1A1) + (0.156*dat$SFRP4)
      androgen_module = (0.634*dat$FAM13C) + (1.079*dat$KLK2) + 
        (0.997*dat$SRD5A2_bounded) + (0.642*dat$AZGP1)
      proliferation_module = dat$TPX2_bounded
      
      risk_score = 0.735*stromal_module - 0.368*cellular_organization_module -
        0.352*androgen_module + 0.095*proliferation_module
      
      return(risk_score)
      
    } else {
      
      stop("The following required Oncotype DX genes: ", 
           paste(setdiff(oncotype_genes, colnames(dat)), collapse = ", "),
           " are not present in colnames of the data")
    }
    
  } else if (base::tolower(test) %in% "decipher") {
    
    if(length(intersect(colnames(dat), decipher_genes)) != 19){
      
      #stop("Only ",
      #     length(intersect(colnames(dat), decipher_genes)),
      #     " out of 19 required genes for Decipher risk score were found")
      warning("The following required Decipher genes: ", 
           paste(setdiff(decipher_genes, colnames(dat)), collapse = ", "),
           " are not present in colnames of the GEX data")
    }
    
    # Overrepresented genes (increase the risk value)
    over <- c("CAMK2N1","EPPK1","IQGAP3","LASP1","NFIB","NUSAP1","PBX1","S1PR4","THBS2","UBE2C","ZWILCH")
    
    # Underrepresented genes (decrease the risk value)
    under <- c("ANO7",
    	"C6orf10","TSBP1", # Aliases
    	"MYBPC1","PCDH7","RABGAP1","TNFRSF19")

    # Intersect between Decipher risk score genes ideally available and current data matrix    
    risk <- dat[,intersect(colnames(dat),c(over, under))]

    # Median centering    
    gene_med <- apply(risk, 2, median)
    risk_centered <- risk - gene_med
    
    # average of the log2 normalized values for the 9 over-expressed targets
    c1 <- apply(risk[,over], MARGIN=1, FUN=function(x) { mean(x, na.rm=TRUE) })
    
    # average of the log2 normalized values for the 9 under-expressed targets
    c2 <- apply(risk[,under], MARGIN=1, FUN=function(x) { mean(x, na.rm=TRUE) })

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
	
	if(base::tolower(test) == "ar"){		
		# Aliases queried using https://www.genecards.org/
		ar_genes <- list(
			"KLK3" = c("KLK3", "PSA", "APS", "KLK2A1"), # Possibly HK3; ambiguous
			"KLK2" = c("KLK2", "HGK-1", "HGK.1", "KLK2A2"), # Possibly HK2; ambiguous
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
	}
}
