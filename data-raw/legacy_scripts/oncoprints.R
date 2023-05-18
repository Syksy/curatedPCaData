###
#
# Oncoprint meta-analyses
# TCGA, Taylor/MSKCC, Barbieri/BROAD, and Ren/EurUrol2017 datasets from GEO/cBioPortal
#
###

# Required libraries, shortens code a lot for not having to use :: and :::
library(MultiAssayExperiment)
library(curatedPCaData)
library(ComplexHeatmap)

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

## TCGA

oncop_tcga <- curatedPCaData:::generate_cbioportal_oncoprint(study_id = "tcga", oncoprintify = TRUE, genes = genes)
# Subset tcga to only the 333 samples
oncop_tcga <- oncop_tcga[, colnames(curatedPCaData::mae_tcga[["gex"]])]

## Taylor et al. / MSKCC

oncop_taylor <- curatedPCaData:::generate_cbioportal_oncoprint(study_id = "taylor", oncoprintify = TRUE, genes = genes)
# Taylor fusions are derived from GEX and will be appended manually from the clinical data
oncop_taylor["ERG", rownames(MultiAssayExperiment::colData(curatedPCaData::mae_taylor))] <-
  ifelse(c("", "Fusion")[MultiAssayExperiment::colData(curatedPCaData::mae_taylor)$ERG_fusion_GEX + 1] == "Fusion",
    # Fusion
    ifelse(oncop_taylor["ERG", rownames(MultiAssayExperiment::colData(curatedPCaData::mae_taylor))] == "",
      # Empty oncoprint slot,
      c("", "Fusion")[MultiAssayExperiment::colData(curatedPCaData::mae_taylor)$ERG_fusion_GEX + 1],
      # Existing genetic aberration at location, appending
      paste(oncop_taylor["ERG", rownames(MultiAssayExperiment::colData(curatedPCaData::mae_taylor))], c("", "Fusion")[MultiAssayExperiment::colData(curatedPCaData::mae_taylor)$ERG_fusion_GEX + 1], sep = ";")
    ),
    # No fusion or NA
    oncop_taylor["ERG", rownames(MultiAssayExperiment::colData(curatedPCaData::mae_taylor))]
  )
# Intersect Taylor et al. to samples with both GEX and CNA based on the sampleMap-function
instances <- as.data.frame(MultiAssayExperiment::sampleMap(curatedPCaData::mae_taylor)) %>%
  dplyr::filter(assay %in% c("cna", "gex")) %>%
  dplyr::select(primary)
# Pick ones with CNA & GEX
oncop_taylor <- oncop_taylor[, names(table(instances)[which(table(instances) >= 2)])]
# Omit genes interesting for fusions which weren't called
oncop_taylor <- oncop_taylor[-which(rownames(oncop_taylor) %in% c("ETV1", "ETV4", "FLI1")), ]

## Barbieri et al./BROAD

oncop_barbieri <- curatedPCaData:::generate_cbioportal_oncoprint(study_id = "barbieri", oncoprintify = TRUE, genes = genes)
# Barbieri fusions for ERG are obtainable from the clinical table
oncop_barbieri["ERG", rownames(MultiAssayExperiment::colData(curatedPCaData::mae_barbieri))] <-
  ifelse(c("", "Fusion")[as.numeric(MultiAssayExperiment::colData(curatedPCaData::mae_barbieri)$TMPRSS2_ERG_FUSION_STATUS %in% c("Positive", "Positive with interstitial deletion")) + 1] == "Fusion",
    # Fusion
    ifelse(oncop_barbieri["ERG", rownames(MultiAssayExperiment::colData(curatedPCaData::mae_barbieri))] == "",
      # Empty oncoprint slot,
      c("", "Fusion")[as.numeric(MultiAssayExperiment::colData(curatedPCaData::mae_barbieri)$TMPRSS2_ERG_FUSION_STATUS %in% c("Positive", "Positive with interstitial deletion")) + 1],
      # Existing genetic aberration at location, appending
      paste(oncop_barbieri["ERG", rownames(MultiAssayExperiment::colData(curatedPCaData::mae_barbieri))], c("", "Fusion")[as.numeric(MultiAssayExperiment::colData(curatedPCaData::mae_barbieri)$TMPRSS2_ERG_FUSION_STATUS %in% c("Positive", "Positive with interstitial deletion")) + 1], sep = ";")
    ),
    # No fusion or NA
    oncop_barbieri["ERG", rownames(MultiAssayExperiment::colData(curatedPCaData::mae_barbieri))]
  )
# Omit genes interesting for fusions which weren't called
oncop_barbieri <- oncop_barbieri[-which(rownames(oncop_barbieri) %in% c("ETV1", "ETV4", "FLI1")), ]
# Subset Barbieri to samples with both CNA (n=109) and MUT (n=112) data
oncop_barbieri <- oncop_barbieri[, intersect(colnames(curatedPCaData::mae_barbieri[["cna"]]), colnames(curatedPCaData::mae_barbieri[["mut"]]))]


## Ren et al. /EurUrol2017

oncop_ren <- curatedPCaData:::generate_cbioportal_oncoprint(study_id = "ren", oncoprintify = TRUE, genes = genes)
# Subset Ren to the mutation called samples (n=63); CNA and GEX available for all (n=65)
#> all(colnames(curatedPCaData::mae_ren[["mut"]]) %in% colnames(curatedPCaData::mae_ren[["gex"]]))
# [1] TRUE
#> all(colnames(curatedPCaData::mae_ren[["mut"]]) %in% colnames(curatedPCaData::mae_ren[["cna"]]))
# [1] TRUE
oncop_ren <- oncop_ren[, colnames(curatedPCaData::mae_ren[["mut"]])]


library(ComplexHeatmap)

## Old functions

# Limit matrix to 'Fusion' events
getFusion <- function(x) {
  tmp <- t(apply(x[genes_fusion, ], MARGIN = 1, FUN = function(z) {
    gsub("HETLOSS|LOW_GAIN|HOMDEL|HIGH_AMP|;|Missense_Mutation|Frame_Shift_Del|Splice_Site|In_Frame_Del|In_Frame_Ins|Frame_Shift_Ins|Nonsense_Mutation", "", z)
  }))
  rownames(tmp) <- paste0("Fusion_", rownames(tmp))
  tmp
}
# Limit matrix to (other) mutation events
getMutations <- function(x) {
  tmp <- t(apply(x[genes_mutations, ], MARGIN = 1, FUN = function(z) {
    gsub("^;", "", gsub("HETLOSS|LOW_GAIN|HOMDEL|HIGH_GAIN|Fusion|", "", z))
  }))
  rownames(tmp) <- paste0("Mutation_", rownames(tmp))
  tmp
}
# Limit matrix to CNA deletions
getDeletions <- function(x) {
  tmp <- t(apply(x[genes_deletions, ], MARGIN = 1, FUN = function(z) {
    gsub("LOW_GAIN|HIGH_AMP|;|Fusion|Missense_Mutation|Frame_Shift_Del|Splice_Site|In_Frame_Del|In_Frame_Ins|Frame_Shift_Ins|Nonsense_Mutation", "", z)
  }))
  rownames(tmp) <- paste0("Deletion_", rownames(tmp))
  tmp
}
# Limit matrix to CNA amplification
getAmplifications <- function(x) {
  tmp <- t(apply(x[genes_amplification, ], MARGIN = 1, FUN = function(z) {
    gsub("HETLOSS|HOMDEL|;|Fusion|Missense_Mutation|Frame_Shift_Del|Splice_Site|In_Frame_Del|In_Frame_Ins|Frame_Shift_Ins|Nonsense_Mutation", "", z)
  }))
  rownames(tmp) <- paste0("Amplification_", rownames(tmp))
  tmp
}

## New functions
# Limit matrix to 'Fusion' events
getFusion <- function(x) {
  tmp <- t(apply(x, MARGIN = 1, FUN = function(z) {
    gsub("HETLOSS|LOW_GAIN|HOMDEL|HIGH_AMP|;|Missense_Mutation|Frame_Shift_Del|Splice_Site|In_Frame_Del|In_Frame_Ins|Frame_Shift_Ins|Nonsense_Mutation|Nonstop_Mutation", "", z)
  }))
  rownames(tmp) <- paste0("Fusion_", rownames(tmp))
  tmp
}
getGain <- function(x) {
  tmp <- t(apply(x, MARGIN = 1, FUN = function(z) {
    gsub("HETLOSS|HOMDEL|Fusion|Missense_Mutation", "", z)
  }))
  rownames(tmp) <- paste0("GainFunction_", rownames(tmp))
  tmp
}
getLoss <- function(x) {
  tmp <- t(apply(x, MARGIN = 1, FUN = function(z) {
    gsub("LOW_GAIN|HIGH_AMP|Fusion", "", z)
  }))
  rownames(tmp) <- paste0("LossFunction_", rownames(tmp))
  tmp
}

GainFunction <- c("MYC", "AKT1", "AR", "PIK3CA")
LossFunction <- c("SPOP", "TP53", "FOXA1", "CDKN1B", "RB1", "MAP3K7", "PTEN", "CDKN1B", "CHD1", "BRCA2", "BRCA1", "ATM", "PARP1")
FusionsFull <- c("ERG", "ETV1", "ETV4", "FLI1")

# TCGA oncoprint with selectivity
tcga <- rbind(
  getFusion(oncop_tcga[FusionsFull, ]),
  getGain(oncop_tcga[intersect(GainFunction, rownames(oncop_tcga)), ]),
  getLoss(oncop_tcga[intersect(LossFunction, rownames(oncop_tcga)), ])
)

# Taylor/MSKCC oncoprint with selectivity
taylor <- rbind(
  Fusion_ERG = gsub("HETLOSS|LOW_GAIN|HOMDEL|HIGH_AMP|;|Missense_Mutation|Frame_Shift_Del|Splice_Site|In_Frame_Del|In_Frame_Ins|Frame_Shift_Ins|Nonsense_Mutation|Nonstop_Mutation", "", oncop_taylor["ERG", ]),
  getGain(oncop_taylor[intersect(GainFunction, rownames(oncop_taylor)), ]),
  getLoss(oncop_taylor[intersect(LossFunction, rownames(oncop_taylor)), ])
)
# Barbieri oncoprint with selectivity
barbieri <- rbind(
  Fusion_ERG = gsub("HETLOSS|LOW_GAIN|HOMDEL|HIGH_AMP|;|Missense_Mutation|Frame_Shift_Del|Splice_Site|In_Frame_Del|In_Frame_Ins|Frame_Shift_Ins|Nonsense_Mutation|Nonstop_Mutation", "", oncop_barbieri["ERG", ]),
  getGain(oncop_barbieri[intersect(GainFunction, rownames(oncop_barbieri)), ]),
  getLoss(oncop_barbieri[intersect(LossFunction, rownames(oncop_barbieri)), ])
)
# Ren oncoprint with selectivity
ren <- rbind(
  Fusion_ERG = gsub("HETLOSS|LOW_GAIN|HOMDEL|HIGH_AMP|;|Missense_Mutation|Frame_Shift_Del|Splice_Site|In_Frame_Del|In_Frame_Ins|Frame_Shift_Ins|Nonsense_Mutation|Nonstop_Mutation", "", oncop_ren["ERG", ]),
  getGain(oncop_ren[intersect(GainFunction, rownames(oncop_ren)), ]),
  getLoss(oncop_ren[intersect(LossFunction, rownames(oncop_ren)), ])
)


# Annotation data
annotations_tcga <- data.frame(
  AR_score = curatedPCaData::mae_tcga[["scores"]]["AR_score", colnames(oncop_tcga)],
  AR_mRNA = curatedPCaData::mae_tcga[["gex"]]["AR", ],
  Gleason = MultiAssayExperiment::colData(curatedPCaData::mae_tcga)[match(colnames(curatedPCaData::mae_tcga[["gex"]]), MultiAssayExperiment::colData(curatedPCaData::mae_tcga)$sample_name), "gleason_grade"],
  Overall_Survival_Status = MultiAssayExperiment::colData(curatedPCaData::mae_tcga)[match(colnames(curatedPCaData::mae_tcga[["gex"]]), MultiAssayExperiment::colData(curatedPCaData::mae_tcga)$sample_name), "overall_survival_status"],
  Days_to_Overall_Survival = MultiAssayExperiment::colData(curatedPCaData::mae_tcga)[match(colnames(curatedPCaData::mae_tcga[["gex"]]), MultiAssayExperiment::colData(curatedPCaData::mae_tcga)$sample_name), "days_to_overall_survival"],
  Disease_Specific_Recurrence_Status = MultiAssayExperiment::colData(curatedPCaData::mae_tcga)[match(colnames(curatedPCaData::mae_tcga[["gex"]]), MultiAssayExperiment::colData(curatedPCaData::mae_tcga)$sample_name), "disease_specific_recurrence_status"],
  Days_to_Disease_Specific_Recurrence = MultiAssayExperiment::colData(curatedPCaData::mae_tcga)[match(colnames(curatedPCaData::mae_tcga[["gex"]]), MultiAssayExperiment::colData(curatedPCaData::mae_tcga)$sample_name), "days_to_disease_specific_recurrence"]
)
annotations_taylor <- data.frame(
  AR_score = mae_taylor[["scores"]]["AR_score", sampleMap(mae_taylor)[sampleMap(mae_taylor)$assay == "gex" & sampleMap(mae_taylor)$primary %in% colnames(oncop_taylor), ]$colname[1:150]],
  AR_mRNA = mae_taylor[["gex"]]["AR", sampleMap(mae_taylor)[sampleMap(mae_taylor)$assay == "gex" & sampleMap(mae_taylor)$primary %in% colnames(oncop_taylor), ]$colname[1:150]],
  Gleason = colData(mae_taylor)[colnames(oncop_taylor), "gleason_grade"],
  Overall_Survival_Status = NA,
  Days_to_Overall_Survival = NA,
  Disease_Specific_Recurrence_Status = colData(mae_taylor)[colnames(oncop_taylor), "disease_specific_recurrence_status"],
  Days_to_Disease_Specific_Recurrence = colData(mae_taylor)[colnames(oncop_taylor), "days_to_disease_specific_recurrence"]
)
rownames(annotations_taylor) <- colnames(oncop_taylor)
annotations_barbieri <- data.frame(
  AR_score = curatedPCaData::mae_barbieri[["scores"]]["AR_score", match(colnames(oncop_barbieri), colnames(curatedPCaData::mae_barbieri[["scores"]]))],
  AR_mRNA = curatedPCaData::mae_barbieri[["gex"]]["AR", match(colnames(oncop_barbieri), colnames(curatedPCaData::mae_barbieri[["gex"]]))],
  Gleason = colData(mae_barbieri)[colnames(oncop_barbieri), "gleason_grade"],
  Overall_Survival_Status = NA,
  Days_to_Overall_Survival = NA,
  Disease_Specific_Recurrence_Status = NA,
  Days_to_Disease_Specific_Recurrence = NA
)
annotations_ren <- data.frame(
  AR_score = curatedPCaData::mae_ren[["scores"]]["AR_score", match(colnames(oncop_ren), colnames(curatedPCaData::mae_ren[["scores"]]))],
  AR_mRNA = curatedPCaData::mae_ren[["gex"]]["AR", match(colnames(oncop_ren), colnames(curatedPCaData::mae_ren[["scores"]]))],
  Gleason = colData(mae_ren)[colnames(oncop_ren), "gleason_grade"],
  Overall_Survival_Status = NA,
  Days_to_Overall_Survival = NA,
  Disease_Specific_Recurrence_Status = NA,
  Days_to_Disease_Specific_Recurrence = NA
)





# ComplexHeatmap colouring options
# Colours for alterations
col <- c(
  "Fusion" = "purple",
  "HOMDEL" = "blue", # Deep deletion
  "HETLOSS" = "lightblue3", # Shallow deletion
  "LOW_GAIN" = "lightpink2", # "Gain"
  "HIGH_AMP" = "red", # "Amplification"
  "Splice_Site" = "greenyellow",
  "Missense_Mutation" = "green4",
  "Nonsense_Mutation" = "green4", # Truncating mutation
  "Frame_Shift_Del" = "orange",
  "Frame_Shift_Ins" = "orange",
  "Nonstop_Mutation" = "green4",
  "In_Frame_Del" = "orange",
  "In_Frame_Ins" = "orange"
)
# Alteration annotation function
alter_fun <- list(
  Fusion = function(x, y, w, h) {
    grid.rect(x, y, w * 0.9, h * 0.7,
      gp = gpar(fill = col["Fusion"], col = NA)
    )
  },
  HOMDEL = function(x, y, w, h) {
    grid.rect(x, y, w * 0.9, h * 0.9,
      gp = gpar(fill = col["HOMDEL"], col = NA)
    )
  },
  HETLOSS = function(x, y, w, h) {
    grid.rect(x, y, w * 0.9, h * 0.9,
      gp = gpar(fill = col["HETLOSS"], col = NA)
    )
  },
  LOW_GAIN = function(x, y, w, h) {
    grid.rect(x, y, w * 0.9, h * 0.9,
      gp = gpar(fill = col["LOW_GAIN"], col = NA)
    )
  },
  HIGH_AMP = function(x, y, w, h) {
    grid.rect(x, y, w * 0.9, h * 0.9,
      gp = gpar(fill = col["HIGH_AMP"], col = NA)
    )
  },
  Splice_Site = function(x, y, w, h) {
    grid.rect(x, y, w * 0.9, h * 0.4,
      gp = gpar(fill = col["Splice_Site"], col = NA)
    )
  },
  Missense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w * 0.9, h * 0.4,
      gp = gpar(fill = col["Missense_Mutation"], col = NA)
    )
  },
  Nonsense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w * 0.9, h * 0.4,
      gp = gpar(fill = col["Nonsense_Mutation"], col = NA)
    )
  },
  Frame_Shift_Del = function(x, y, w, h) {
    grid.rect(x, y, w * 0.9, h * 0.4,
      gp = gpar(fill = col["Frame_Shift_Del"], col = NA)
    )
  },
  Frame_Shift_Ins = function(x, y, w, h) {
    grid.rect(x, y, w * 0.9, h * 0.4,
      gp = gpar(fill = col["Frame_Shift_Ins"], col = NA)
    )
  },
  Nonstop_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w * 0.9, h * 0.4,
      gp = gpar(fill = col["Nonstop_Mutation"], col = NA)
    )
  },
  In_Frame_Del = function(x, y, w, h) {
    grid.rect(x, y, w * 0.9, h * 0.4,
      gp = gpar(fill = col["In_Frame_Del"], col = NA)
    )
  },
  In_Frame_Ins = function(x, y, w, h) {
    grid.rect(x, y, w * 0.9, h * 0.4,
      gp = gpar(fill = col["In_Frame_Ins"], col = NA)
    )
  }
)
# Custom heatmap_legend_param
legpar <- list(
  title = "Alterations",
  at = c("Fusion", "HOMDEL", "HETLOSS", "LOW_GAIN", "HIGH_AMP", "Splice_Site", "Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Nonstop_Mutation", "In_Frame_Del", "In_Frame_Ins"),
  labels = c("Fusion", "Deep Deletion", "Shallow Deletion", "Gain", "Amplification", "Splice Site", "SNV", "SNV", "Indel", "Indel", "SNV", "Indel", "Indel")
)

# Function for counting alterations for a given gene, takes into account
order_oncop <- function(oncop, genes = c("ERG", "MAP3K7", "CHD1", "PTEN", "TP53", "SPOP", "RB1")) {
  # Intersect with available genes
  genes <- intersect(genes, rownames(oncop))
  # Binarization function for oncoprint matrix; alteration (=1) or not (=0)
  binarize <- function(mat) {
    bin <- matrix(as.numeric(!mat == ""), nrow = nrow(mat), ncol = ncol(mat))
    dimnames(bin) <- dimnames(mat)
    bin
  }
  # Genes in which order we'll perform initial sorting
  oncomat <- binarize(oncop)

  # rest_genes <- rownames(oncomat)[!rownames(oncomat) %in% genes]
  # burden <- apply(oncomat[rest_genes,], MARGIN=2, FUN=sum)
  # Rest of gene alterations, will serve as a tiebreaker when sorting the rest
  # oncop <- oncop[,order(burden, decreasing=TRUE)]
  # Re-binarize due to re-organizing
  # Genes in which order we'll perform initial sorting
  oncomat <- binarize(oncop)

  # Order in reverse order so the right genes get assigned left-wards
  ord <- rev(do.call("order", lapply(genes, FUN = function(g) {
    oncomat[g, ]
  })))
  # for(g in rev(genes)){
  # 	oncop <- oncop[,order(oncomat[g,], decreasing=TRUE)]
  # }
  # Return column-ordered oncoprint
  oncop[, ord]
}

# Try to order genes (depending on the intersect, may vary by dataset)
gene_ordering <- c(
  # Fusions
  "ERG", "ETV1", "ETV4", "FLI1",
  # Deletions
  "MAP3K7", "CHD1", "SPOP", "PTEN", "TP53", "RB1", "FOXA1", "CDKN1B", "PARP1", "ATM", "BRCA1", "BRCA2",
  # Amps/gains
  "MYC", "PIK3CA", "AKT1", "AR"
)
# Tidy up prefixes
remove_prefix <- function(x) {
  gsub("Fusion_|GainFunction_|LossFunction_", "", x)
}
rownames(tcga) <- remove_prefix(rownames(tcga))
rownames(taylor) <- remove_prefix(rownames(taylor))
rownames(barbieri) <- remove_prefix(rownames(barbieri))
rownames(ren) <- remove_prefix(rownames(ren))
# Reorder matrices
tcga <- tcga[intersect(gene_ordering, rownames(tcga)), ]
taylor <- taylor[intersect(gene_ordering, rownames(taylor)), ]
barbieri <- barbieri[intersect(gene_ordering, rownames(barbieri)), ]
ren <- ren[intersect(gene_ordering, rownames(ren)), ]

pdf("Oncoprint_TCGA_draft_May27.pdf", width = 12, height = 4)
# ComplexHeatmap OncoPrint for TCGA
ht_tcga <- ComplexHeatmap::oncoPrint(
  # Ordering to the desired priorization of genes
  # tcga,
  order_oncop(tcga),
  column_title = "TCGA",
  row_title = NULL,
  alter_fun = alter_fun,
  row_order = 1:nrow(tcga), # Default gene ordering after taking intersection; do not reorder based on alteration percentages
  column_order = 1:ncol(tcga),
  top_annotation = NULL, right_annotation = NULL, # Omit default annotations
  heatmap_legend_param = legpar,
  col = col
)
ComplexHeatmap::draw(ht_tcga)
dev.off()

# Samples in taylor et all with no 'omics overlap, omit
taylor <- taylor[, -which(apply(taylor, MARGIN = 2, FUN = function(z) {
  all(is.na(z))
}))]
# SPOP-gene spuriously shows too low alteration percentage, as there were a lack of samples with SNV-like mutations called
taylor <- taylor[-which(rownames(taylor) == "SPOP"), ]
pdf("Oncoprint_Taylor_draft_May27.pdf", width = 12, height = 4)
# ComplexHeatmap OncoPrint for Taylor et al. / MSKCC
ht_taylor <- ComplexHeatmap::oncoPrint(
  # taylor,
  order_oncop(taylor),
  column_title = "MSKCC",
  row_title = NULL,
  alter_fun = alter_fun,
  row_order = 1:nrow(taylor), # Default gene ordering after taking intersection; do not reorder based on alteration percentages
  column_order = 1:ncol(taylor),
  top_annotation = NULL, right_annotation = NULL, # Omit default annotations
  heatmap_legend_param = legpar,
  col = col
)
ComplexHeatmap::draw(ht_taylor)
dev.off()

# AR expression not properly reported in barbieri
barbieri <- barbieri[-grep("AR", rownames(barbieri)), ]
pdf("Oncoprint_Barbieri_draft_May27.pdf", width = 12, height = 4)
# ComplexHeatmap OncoPrint for Barbieri et al. / BROAD
ht_barbieri <- ComplexHeatmap::oncoPrint(
  # barbieri,
  order_oncop(barbieri),
  column_title = "BROAD",
  row_title = NULL,
  alter_fun = alter_fun,
  row_order = 1:nrow(barbieri), # Default gene ordering after taking intersection; do not reorder based on alteration percentages
  column_order = 1:ncol(barbieri),
  top_annotation = NULL, right_annotation = NULL, # Omit default annotations
  heatmap_legend_param = legpar,
  col = col
)
ComplexHeatmap::draw(ht_barbieri)
dev.off()

## Omitting Ren / SMMU dataset due to its moderately small size

## ERG fusions not properly reported in ren
# ren <- ren[-grep("ERG", rownames(ren)),]
# pdf("Oncoprint_Ren_draft_May27.pdf", width=12, height=4)
## ComplexHeatmap OncoPrint for Ren et al. / EurUrol2017
# ht_ren <- ComplexHeatmap::oncoPrint(
# 	ren,
# 	column_title = "SMMU",
# 	row_title = NULL,
# 	alter_fun = alter_fun,
# 	row_order = 1:nrow(ren), # Default gene ordering after taking intersection; do not reorder based on alteration percentages
# 	top_annotation = NULL, right_annotation = NULL, # Omit default annotations
# 	heatmap_legend_param = legpar,
# 	col=col
# )
# ComplexHeatmap::draw(ht_ren)
# dev.off()

# Write annotations & raw oncoprints as TSVs
# set to FALSE to comment out for now so whole .R-file script can be run just for heatmaps
if (FALSE) {
  write.table(cbind(annotations_tcga, t(oncop_tcga)), file = "data_tcga.tsv", sep = "\t", quote = FALSE, dec = ",")
  write.table(cbind(annotations_taylor, t(oncop_taylor)), file = "data_taylor.tsv", sep = "\t", quote = FALSE, dec = ",")
  write.table(cbind(annotations_barbieri, t(oncop_barbieri)), file = "data_barbieri.tsv", sep = "\t", quote = FALSE, dec = ",")
  write.table(cbind(annotations_ren, t(oncop_ren)), file = "data_ren.tsv", sep = "\t", quote = FALSE, dec = ",")
}

## N counts for reporting
#> dim(tcga)
# [1]  20 333
#> dim(taylor)
# [1]  17 122
#> dim(barbieri)
# [1]  15 109
#> dim(ren)
# [1] 16 63
