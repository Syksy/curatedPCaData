# Bioconductor;
# Warning! Old packages may cause some dependencies to fail, while their updating may fail if they're already in use for the session
# Safest is to update all Bioconductor packages before analyses / running pipelines

# BioConductor 'annotate' package for all sorts of conversions and genetic location info, etc
# Annotation for microarrays
# For Entrez <-> Hugo Gene Symbol mapping
# http://bioconductor.org/packages/release/bioc/html/annotate.html

# Genome wide annotation for Human
# For Entrez <-> Hugo Gene Symbol mapping (database)
# https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html

# See if legacy aliases are supported
# grep("alias", biomaRt::listAttributes(biomaRt::useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl"))[,"name"], value=TRUE)

# Example for extracting the affy specific probes as suggested in https://www.biostars.org/p/402677/
# mart <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL")
# mart <- biomaRt::useDataset(grep("hsapiens", listDatasets(mart)[,1], value=TRUE), mart)
# listAttributes(mart)[grep("affy", listAttributes(mart)[,1]),]
# listAttributes(mart)[grep("agilent", listAttributes(mart)[,1]),]
#> listAttributes(mart)[grep("affy", listAttributes(mart)[,1]),]
#                     name                 description         page
#107          affy_hc_g110          AFFY HC G110 probe feature_page
#108         affy_hg_focus         AFFY HG Focus probe feature_page
#109         affy_hg_u133a         AFFY HG U133A probe feature_page
#110       affy_hg_u133a_2       AFFY HG U133A 2 probe feature_page
#111         affy_hg_u133b         AFFY HG U133B probe feature_page
#112   affy_hg_u133_plus_2   AFFY HG U133 Plus 2 probe feature_page
#113          affy_hg_u95a          AFFY HG U95A probe feature_page
#114        affy_hg_u95av2        AFFY HG U95Av2 probe feature_page
#115          affy_hg_u95b          AFFY HG U95B probe feature_page
#116          affy_hg_u95c          AFFY HG U95C probe feature_page
#117          affy_hg_u95d          AFFY HG U95D probe feature_page
#118          affy_hg_u95e          AFFY HG U95E probe feature_page
#119          affy_hta_2_0          AFFY HTA 2 0 probe feature_page
#120   affy_huex_1_0_st_v2   AFFY HuEx 1 0 st v2 probe feature_page
#121         affy_hugenefl         AFFY HuGeneFL probe feature_page
#122 affy_hugene_1_0_st_v1 AFFY HuGene 1 0 st v1 probe feature_page
#123 affy_hugene_2_0_st_v1 AFFY HuGene 2 0 st v1 probe feature_page
#124        affy_primeview        AFFY PrimeView probe feature_page
#125         affy_u133_x3p         AFFY U133 X3P probe feature_page
#> listAttributes(mart)[grep("agilent", listAttributes(mart)[,1]),]
#                                name                            description         page
#126                  agilent_cgh_44b                  AGILENT CGH 44b probe feature_page
#127                 agilent_gpl26966                 AGILENT GPL26966 probe feature_page
#128                  agilent_gpl6848                  AGILENT GPL6848 probe feature_page
#129    agilent_sureprint_g3_ge_8x60k    AGILENT SurePrint G3 GE 8x60k probe feature_page
#130 agilent_sureprint_g3_ge_8x60k_v2 AGILENT SurePrint G3 GE 8x60k v2 probe feature_page
#131              agilent_wholegenome              AGILENT WholeGenome probe feature_page
#132     agilent_wholegenome_4x44k_v1     AGILENT WholeGenome 4x44k v1 probe feature_page
#133     agilent_wholegenome_4x44k_v2     AGILENT WholeGenome 4x44k v2 probe feature_page  
  
mart <- biomaRt::useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl"  
# Fetch gene names for various annotations
# Gene type annotations
curatedPCaData_genes <- biomaRt::getBM(
	attributes = 
		c(
			# Hugo
			'hgnc_symbol',
			# ENSEMBL
			'ensembl_gene_id', 'ensembl_transcript_id',
			# RefSeq
			'refseq_mrna',
			# Chromosomal information
			'chromosome_name','start_position','end_position',
			# Description
			'description'
		),
	mart = mart)
)
# Annotations based on specific array platforms
curatedPCaData_genes_affy_hg_u133a <- biomaRt::getBM(
	attributes = 
		c(
			# Hugo
			'hgnc_symbol',
			# Array specific probe identifiers needed by processed studies
			"affy_hg_u133a"	# Used by: Sun et al., ...
		),
	mart = mart)
)
curatedPCaData_genes_affy_hg_u133a_2 <- biomaRt::getBM(
	attributes = 
		c(
			# Hugo
			'hgnc_symbol',
			# Array specific probe identifiers needed by processed studies
			"affy_hg_u133a_2"	# Used by: Wallace et al., IGC... 
		),
	mart = mart)
)
# Timing out
curatedPCaData_genes_affy_huex_1_0_st_v2 <- biomaRt::getBM(
	attributes = 
		c(
			# Hugo
			'hgnc_symbol',
			# Array specific probe identifiers needed by processed studies
			"affy_huex_1_0_st_v2"	# Used by: Taylor et al., ...
		),
	mart = mart)
)

# Sort the generic list first using chromosomes and then bp locations
curatedPCaData_genes <- curatedPCaData_genes[order(curatedPCaData_genes$chromosome_name, curatedPCaData_genes$start_position, curatedPCaData_genes$end_position),]
# Connect to annotation database and extract a list of gene name synonyms, aliases, legacy names
db.con <- org.Hs.eg.db::org.Hs.eg_dbconn()
# SQL query 
aliases <- DBI::dbGetQuery(db.con, "SELECT * FROM alias, gene_info WHERE alias._id == gene_info._id;")
curatedPCaData_genes[,"Aliases"] <- unlist(lapply(curatedPCaData_genes$hgnc_symbol, FUN=function(g){
	w <- which(aliases[,"symbol"] == g)
	if(length(w)>0){
		paste0(aliases[w,"alias_symbol"], collapse=";")
	}else{
		""
	}
}))
# Omit row names (wrong order indices)
rownames(curatedPCaData_genes) <- NULL
# Save gene information extraction date as an attribute 'date'
attr(curatedPCaData_genes, 'date') <- Sys.time()
# Save the data frame of various gene annotations for package's internal use
usethis::use_data(
	# Unquoted names for objects to save
	curatedPCaData_genes, curatedPCaData_genes_affy_hg_u133a, curatedPCaData_genes_affy_hg_u133a_2,
	## Timeout issues: curatedPCaData_genes_affy_huex_1_0_st_v2
	# Saving parameters
	internal = TRUE, overwrite = TRUE)