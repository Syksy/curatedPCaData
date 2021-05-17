###
#
# Oncoprint meta-analyses
# TCGA, Taylor/MSKCC, Barbieri/BROAD, and Ren/EurUrol2017 datasets from GEO/cBioPortal
#
###

genes_fusion <- c("ERG", "ETV1", "ETV4", "FLI1")
genes_mutations <- c("SPOP", "TP53", "FOXA1", "AKT1", "CDKN1B", "PIK3CA", "RB1", "BRCA1", "BRCA2")
genes_deletions <- c("MAP3K7", "CHD1", "RB1", "TP53", "PTEN", "BRCA1", "BRCA2", "PARP1", "ATM", "CDKN1B")
genes_amplification <- c("MYC", "AR", "AKT1")

# List for all genes of interest together
genes <- unique(c(genes_fusion, genes_mutations, genes_deletions, genes_amplification))

## Meta-data to plot
# - AR-score
# - AR mRNA expression (z-score)
# - Cohort
# - Gleason score
# - Recurrence/Treatments/Overall Survival

oncop_tcga <- curatedPCaData:::generate_cbioportal_oncoprint(study_id="tcga", oncoprintify=TRUE, genes=genes)
oncop_taylor <- curatedPCaData:::generate_cbioportal_oncoprint(study_id="taylor", oncoprintify=TRUE, genes=genes)
oncop_barbieri <- curatedPCaData:::generate_cbioportal_oncoprint(study_id="barbieri", oncoprintify=TRUE, genes=genes)
oncop_ren <- curatedPCaData:::generate_cbioportal_oncoprint(study_id="ren", oncoprintify=TRUE, genes=genes)


library(ComplexHeatmap)

# Limit matrix to 'Fusion' events
getFusion <- function(x){
	tmp <- t(apply(x[genes_fusion,], MARGIN=1, FUN=function(z) { gsub("HETLOSS|LOW_GAIN|HOMDEL|HIGH_AMP|;|Missense_Mutation|Frame_Shift_Del|Splice_Site|In_Frame_Del|In_Frame_Ins|Frame_Shift_Ins|Nonsense_Mutation", "", z) }))
	rownames(tmp) <- paste0("Fusion_", rownames(tmp))
	tmp
}
# Limit matrix to (other) mutation events
getMutations <- function(x){
	tmp <- t(apply(x[genes_mutations,], MARGIN=1, FUN=function(z) { gsub("^;", "", gsub("HETLOSS|LOW_GAIN|HOMDEL|HIGH_GAIN|Fusion|", "", z)) }))
	rownames(tmp) <- paste0("Mutation_", rownames(tmp))
	tmp
}
# Limit matrix to CNA deletions
getDeletions <- function(x){
	tmp <- t(apply(x[genes_deletions,], MARGIN=1, FUN=function(z) { gsub("LOW_GAIN|HIGH_AMP|;|Fusion|Missense_Mutation|Frame_Shift_Del|Splice_Site|In_Frame_Del|In_Frame_Ins|Frame_Shift_Ins|Nonsense_Mutation", "", z) }))
	rownames(tmp) <- paste0("Deletion_", rownames(tmp))
	tmp
}
# Limit matrix to CNA amplification
getAmplifications <- function(x){
	tmp <- t(apply(x[genes_amplification,], MARGIN=1, FUN=function(z) { gsub("HETLOSS|HOMDEL|;|Fusion|Missense_Mutation|Frame_Shift_Del|Splice_Site|In_Frame_Del|In_Frame_Ins|Frame_Shift_Ins|Nonsense_Mutation", "", z) }))
	rownames(tmp) <- paste0("Amplification_", rownames(tmp))
	tmp
}

# TCGA oncoprint with selectivity
tcga <- rbind(
	getFusion(oncop_tcga[genes_fusion,colnames(curatedPCaData::mae_tcga[["cna"]])]),
	getMutations(oncop_tcga[genes_mutations,colnames(curatedPCaData::mae_tcga[["cna"]])]),
	getDeletions(oncop_tcga[genes_deletions,colnames(curatedPCaData::mae_tcga[["cna"]])]),
	getAmplifications(oncop_tcga[genes_amplification,colnames(curatedPCaData::mae_tcga[["cna"]])])
)
# Taylor/MSKCC oncoprint with selectivity
taylor <- rbind(
	getFusion(oncop_taylor[genes_fusion,]),
	getMutations(oncop_taylor[genes_mutations,]),
	getDeletions(oncop_taylor[genes_deletions,]),
	getAmplifications(oncop_taylor[genes_amplification,])
)
# Barbieri oncoprint with selectivity
barbieri <- rbind(
	getFusion(oncop_barbieri[genes_fusion,]),
	getMutations(oncop_barbieri[genes_mutations,]),
	getDeletions(oncop_barbieri[genes_deletions,]),
	getAmplifications(oncop_barbieri[genes_amplification,])
)
# Ren oncoprint with selectivity
ren <- rbind(
	getFusion(oncop_ren[genes_fusion,]),
	getMutations(oncop_ren[genes_mutations,]),
	getDeletions(oncop_ren[genes_deletions,]),
	getAmplifications(oncop_ren[genes_amplification,])
)





# Alteration annotation function
alter_fun = list(
	Fusion = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
		gp = gpar(fill = col["Fusion"], col = NA)),
	HOMDEL = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, 
		gp = gpar(fill = col["HOMDEL"], col = NA)),
	HETLOSS = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, 
		gp = gpar(fill = col["HETLOSS"], col = NA)),
	LOW_GAIN = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, 
		gp = gpar(fill = col["LOW_GAIN"], col = NA)),
	HIGH_AMP = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, 
		gp = gpar(fill = col["HIGH_AMP"], col = NA)),
	Splice_Site = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, 
		gp = gpar(fill = col["Splice_Site"], col = NA)),
	Missense_Mutation = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, 
		gp = gpar(fill = col["Missense_Mutation"], col = NA)),
	Nonsense_Mutation = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, 
		gp = gpar(fill = col["Nonsense_Mutation"], col = NA)),
	Frame_Shift_Del = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, 
		gp = gpar(fill = col["Frame_Shift_Del"], col = NA)),
	Frame_Shift_Ins = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, 
		gp = gpar(fill = col["Frame_Shift_Ins"], col = NA)),
	Nonstop_Mutation = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, 
		gp = gpar(fill = col["Nonstop_Mutation"], col = NA)),
	In_Frame_Del = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, 
		gp = gpar(fill = col["In_Frame_Del"], col = NA)),
	In_Frame_Ins = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, 
		gp = gpar(fill = col["In_Frame_Ins"], col = NA))
)
# Colours for alterations
col = c(
	"Fusion" = rainbow(13)[1],
	"HOMDEL" = rainbow(13)[2],
	"HETLOSS" = rainbow(13)[3],
	"LOW_GAIN" = rainbow(13)[4],
	"HIGH_AMP" = rainbow(13)[5],
	"Splice_Site" = rainbow(13)[6],
	"Missense_Mutation" = rainbow(13)[7],
	"Nonsense_Mutation" = rainbow(13)[8],
	"Frame_Shift_Del" = rainbow(13)[9],
	"Frame_Shift_Ins" = rainbow(13)[10],
	"Nonstop_Mutation" = rainbow(13)[11],
	"In_Frame_Del" = rainbow(13)[12],
	"In_Frame_Ins" = rainbow(13)[13]
)

# Rows splitting in each heatmap
rowsplits <- c(
	rep("1) Fusions", times=length(genes_fusion)),
	rep("2) Mutations", times=length(genes_mutations)),
	rep("3) Deletions", times=length(genes_deletions)),
	rep("4) Amplifications", times=length(genes_amplification))
)

# ComplexHeatmap OncoPrint for TCGA	
chop_tcga <- ComplexHeatmap::oncoPrint(
	tcga, 
	column_title = "TCGA",
	alter_fun = alter_fun, 
	col=col,
	row_split = rowsplits,
	show_row_names = FALSE,
	right_annotation = NULL,	
	# Top annotations
	top_annotation = ComplexHeatmap::HeatmapAnnotation(
		"AR score" = curatedPCaData::mae_tcga[["scores"]]["AR_score",colnames(tcga)],
		"AR mRNA" = curatedPCaData::mae_tcga[["gex"]]["AR",colnames(tcga)],
		"Gleason grade" = MultiAssayExperiment::colData(curatedPCaData::mae_tcga)[match(colnames(tcga), MultiAssayExperiment::colData(curatedPCaData::mae_tcga)$sample_name),"gleason_grade"],
		col = list(
			"Gleason grade" = c("5" = "lightgreen", "6" = "green", "7" = "cyan", "8" = "blue", "9" = "yellow", "10" = "orange")
		)
	)
)
# ComplexHeatmap OncoPrint for Taylor et al./MSKCC	
chop_taylor <- ComplexHeatmap::oncoPrint(
	taylor, 
	column_title = "Taylor/MSKCC",	
	alter_fun = alter_fun, 
	col=col,
	row_split = rowsplits,
	show_row_names = FALSE,
	right_annotation = NULL,	
	# Top annotations
	top_annotation = ComplexHeatmap::HeatmapAnnotation(
		"AR score" = curatedPCaData::mae_taylor[["scores"]]["AR_score",match(colnames(taylor), colnames(curatedPCaData::mae_taylor[["scores"]]))],
		"AR mRNA" = curatedPCaData::mae_taylor[["gex"]]["AR",match(colnames(taylor), colnames(curatedPCaData::mae_taylor[["gex"]]))],
		"Gleason grade" = MultiAssayExperiment::colData(curatedPCaData::mae_taylor)[match(colnames(taylor), MultiAssayExperiment::colData(curatedPCaData::mae_taylor)$patient_id),"gleason_grade"],
		col = list(
			"Gleason grade" = c("5" = "lightgreen", "6" = "green", "7" = "cyan", "8" = "blue", "9" = "yellow", "10" = "orange")
		)
	)
)
# ComplexHeatmap OncoPrint for Barbieri et al.
chop_barbieri <- ComplexHeatmap::oncoPrint(
	barbieri, 
	column_title = "Barbieri/BROAD",	
	alter_fun = alter_fun, 
	col=col,
	row_split = rowsplits,
	show_row_names = FALSE,
	right_annotation = NULL,	
	# Top annotations
	top_annotation = ComplexHeatmap::HeatmapAnnotation(
		"AR score" = curatedPCaData::mae_barbieri[["scores"]]["AR_score",match(colnames(barbieri), colnames(curatedPCaData::mae_barbieri[["scores"]]))],
		"AR mRNA" = curatedPCaData::mae_barbieri[["gex"]]["AR",match(colnames(barbieri), colnames(curatedPCaData::mae_barbieri[["gex"]]))],
		"Gleason grade" = MultiAssayExperiment::colData(curatedPCaData::mae_barbieri)[match(colnames(barbieri), MultiAssayExperiment::colData(curatedPCaData::mae_barbieri)$sample_name),"gleason_grade"],
		col = list(
			"Gleason grade" = c("3+3" = "green", "4+3" = "cyan", "3+4" = "cyan", "3+5" = "blue", "4+4" = "blue", "4+5" = "yellow", "5+4" = "yellow", "5+5" = "orange")
		)
	)
)
# ComplexHeatmap OncoPrint for Ren et al.
chop_ren <- ComplexHeatmap::oncoPrint(
	ren, 
	column_title = "Ren/EurUrol2017",	
	alter_fun = alter_fun, 
	col=col,
	row_split = rowsplits,
	show_row_names = TRUE,
	right_annotation = NULL,	
	# Top annotations
	top_annotation = ComplexHeatmap::HeatmapAnnotation(
		"AR score" = curatedPCaData::mae_ren[["scores"]]["AR_score",colnames(ren)],
		"AR mRNA" = curatedPCaData::mae_ren[["gex"]]["AR",colnames(ren)],
		"Gleason grade" = MultiAssayExperiment::colData(curatedPCaData::mae_ren)[match(colnames(ren), MultiAssayExperiment::colData(curatedPCaData::mae_ren)$sample_name),"gleason_grade"],
		col = list(
			"Gleason grade" = c("3+3" = "green", "4+3" = "cyan", "3+4" = "cyan", "3+5" = "blue", "4+4" = "blue", "4+5" = "yellow", "5+4" = "yellow", "5+5" = "orange")
		)
	)
)

ComplexHeatmap::draw(chop_tcga + chop_taylor + chop_barbieri + chop_ren)
