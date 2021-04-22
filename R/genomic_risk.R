#' Genomic risk scores
#'
#' Prolaris, Oncotype DX, Decipher genomic panels for prostate cancer risk
#'
#' @details https://bjui-journals.onlinelibrary.wiley.com/doi/10.1111/bju.14452 https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-14-690 https://www.nature.com/articles/s41391-019-0167-9
#'
genomic_risk <- function(mae, 
                         object = "gex",
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
  
  
  dat <- mae@ExperimentList[[object]]
  dat <- t(dat) 
  dat <- as.data.frame(dat)
  
  if(test == "Prolaris"){
    
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
    
  } else if (test %in% c("Oncotype DX", "OncotypeDX")){
    
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
    
  } else if (test %in% "Decipher") {
    
    if(length(intersect(colnames(dat), decipher_genes)) != 19){
      
      #stop("Only ",
      #     length(intersect(colnames(dat), decipher_genes)),
      #     " out of 19 required genes for Decipher risk score were found")
      stop("The following required Decipher genes: ", 
           paste(setdiff(decipher_genes, colnames(dat)), collapse = ", "),
           " are not present in colnames of the data")
    }
    
    risk <- dat[,colnames(dat) %in% decipher_genes]
    gene_med <- apply(risk, 2, median)
    risk_centered <- risk - gene_med
    
    # cluster1 <- hclust(dist(t(risk)), method="centroid")
    
    over <- c("CAMK2N1","EPPK1","IQGAP3","LASP1","NFIB","NUSAP1","PBX1","S1PR4","THBS2","UBE2C",
              "ZWILCH")
    
    under <- c("ANO7","C6orf10","MYBPC1","PCDH7","RABGAP1","TNFRSF19")
    
    # average of the log2 normalized values for the 9 over-expressed targets
    c1 <- apply(risk[,c("CAMK2N1","EPPK1","IQGAP3","LASP1","NFIB","NUSAP1","PBX1",
                        "S1PR4","THBS2","UBE2C","ZWILCH")], 1, mean)
    
    # average of the log2 normalized values for the 9 under-expressed targets
    c2 <- apply(risk[,c("ANO7","C6orf10","MYBPC1","PCDH7","RABGAP1","TNFRSF19")], 1, mean)
    
    risk_score <- c1-c2
    
    return(risk_score)
    
  } else {
  	stop(paste("Invalid genomic risk score name:", test))
  }
}