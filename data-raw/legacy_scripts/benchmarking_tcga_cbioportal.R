tcga_firehouse <- rio::import("/Users/varsha/Downloads/prad_tcga/data_RNA_Seq_v2_expression_median.txt")
tcga_firehouse <- tcga_firehouse[, -2]
names(tcga_firehouse) <- gsub(x = names(tcga_firehouse), pattern = "\\-", replacement = ".")
tcga_firehouse[, 2:499] <- as.numeric(unlist(tcga_firehouse[, 2:499]))
tcga_firehouse2 <- tcga_firehouse[-c(3, 3056, 3185, 4619, 8940, 9775, 11582, 11751, 11891, 12055, 12899, 13597, 14474, 15958, 18906), ]
tcga_firehouse2 <- as.data.frame(tcga_firehouse2)
tcga_firehouse2 <- tcga_firehouse2[complete.cases(tcga_firehouse2[, 1]), ]
rownames(tcga_firehouse2) <- tcga_firehouse2[, 1]
tcga_firehouse2 <- tcga_firehouse2[, -1]

tcga_pancancer <- rio::import("/Users/varsha/Downloads/prad_tcga_pan_can_atlas_2018/data_RNA_Seq_v2_expression_median.txt")
tcga_pancancer <- tcga_pancancer[, -2]
names(tcga_pancancer) <- gsub(x = names(tcga_pancancer), pattern = "\\-", replacement = ".")
tcga_pancancer <- tcga_pancancer[!grepl("-", tcga_pancancer$Hugo_Symbol), ]
tcga_pancancer[, 2:494] <- as.numeric(unlist(tcga_pancancer[, 2:494]))
tcga_pancancer <- as.data.frame(tcga_pancancer)
tcga_pancancer <- tcga_pancancer[!(is.na(tcga_pancancer$Hugo_Symbol) | tcga_pancancer$Hugo_Symbol == ""), ]

# tcga_pancancer<- tcga_pancancer[-c(1,2,5,6,7,10,14,19:22,25,27),]
a <- duplicated(tcga_pancancer$Hugo_Symbol)
a <- as.data.frame(a)
tcga_pancancer <- tcga_pancancer[-c(9652, 11647, 12640, 14211, 17976), ]

O <- apply(n, 1, function(x) any(x %in% "TRUE"))
O <- as.data.frame(O)
write.csv(O, "/Users/varsha/Desktop/O.csv")
# tcga_pancancer2<- tcga_pancancer[-c(11498:11512,11539),]
# tcga_pancancer2<-tcga_pancancer2[-grep("-", tcga_pancancer2$Hugo_Symbol),]
#
# tcga_pancancer2<- tcga_pancancer2[-11539,]
# #tcga_pancancer2$Hugo_Symbol<-gsub("-","",as.character(tcga_pancancer2$Hugo_Symbol))
# b <- duplicated(tcga_pancancer2$Hugo_Symbol)
# b <- as.data.frame(b)
# tcga_pancancer2<- tcga_pancancer2[-12883,]

rownames(tcga_pancancer) <- tcga_pancancer$Hugo_Symbol
tcga_pancancer <- tcga_pancancer[, -1]


# a<- duplicated(tcga_pancancer[,1])
# a<- as.data.frame(a)

tcga_cell_2015 <- curatedPCaData::mae_tcga[["gex"]]


library(immunedeconv)

epic_tcga_firehouse <- immunedeconv::deconvolute(tcga_firehouse2, method = "epic")
quantiseq_tcga_firehouse <- immunedeconv::deconvolute(tcga_firehouse2, method = "quantiseq")
xcell_tcga_firehouse <- immunedeconv::deconvolute(tcga_firehouse2, method = "xcell")
mcp_tcga_firehouse <- immunedeconv::deconvolute(tcga_firehouse2, method = "mcp_counter")

epic_tcga_pancancer <- immunedeconv::deconvolute(tcga_pancancer, method = "epic")
quantiseq_tcga_pancancer <- immunedeconv::deconvolute(tcga_pancancer2, method = "quantiseq")
xcell_tcga_pancancer <- immunedeconv::deconvolute(tcga_pancancer, method = "xcell")
mcp_tcga_pancancer <- immunedeconv::deconvolute(tcga_pancancer, method = "mcp_counter")

# FORMATTING EPIC DATA

epic_tcga_firehouse_t <- t(epic_tcga_firehouse)
colnames(epic_tcga_firehouse_t) <- epic_tcga_firehouse_t[1, ]
epic_tcga_firehouse_t <- epic_tcga_firehouse_t[-1, ]
epic_tcga_firehouse_t <- as.data.frame(epic_tcga_firehouse_t)
epic_tcga_firehouse_t[, 9] <- rownames(epic_tcga_firehouse_t)

epic_tcga_pancancer_t <- t(epic_tcga_pancancer)
colnames(epic_tcga_pancancer_t) <- epic_tcga_pancancer_t[1, ]
epic_tcga_pancancer_t <- epic_tcga_pancancer_t[-1, ]
epic_tcga_pancancer_t <- as.data.frame(epic_tcga_pancancer_t)
epic_tcga_pancancer_t[, 9] <- rownames(epic_tcga_pancancer_t)

######################################################################################
# EPIC - FIREHOUSE VS PANCANCER
######################################################################################

merged <- merge(epic_tcga_firehouse_t, epic_tcga_pancancer_t, by.x = "V9", by.y = "V9")
pdf("data-raw/benchmarking/tcga_firehouse_vs_pancancer_EPIC/Bcells.pdf")
plot(merged$`B cell.x`, merged$`B cell.y`, xlab = "epic_tcga_firehouse_Bcells", ylab = "epic_tcga_pancancer_Bcells")
dev.off()
pdf("data-raw/benchmarking/tcga_firehouse_vs_pancancer_EPIC/Cancer_associated_fibroblasts.pdf")
plot(merged$`Cancer associated fibroblast.x`, merged$`Cancer associated fibroblast.y`, xlab = "epic_tcga_firehouse_CAF", ylab = "epic_tcga_pancancer_CAF")
dev.off()
pdf("data-raw/benchmarking/tcga_firehouse_vs_pancancer_EPIC/T_cell_CD4+.pdf")
plot(merged$`T cell CD4+.x`, merged$`T cell CD4+.y`, xlab = "epic_tcga_firehouse_T_cell_CD4+", ylab = "epic_tcga_pancancer_T_cell_CD4+")
dev.off()
pdf("data-raw/benchmarking/tcga_firehouse_vs_pancancer_EPIC/T_cell_CD8+.pdf")
plot(merged$`T cell CD8+.x`, merged$`T cell CD8+.y`, xlab = "epic_tcga_firehouse_T_cell_CD8+", ylab = "epic_tcga_pancancer_T_cell_CD8+")
dev.off()
pdf("data-raw/benchmarking/tcga_firehouse_vs_pancancer_EPIC/Endothelial_cells.pdf")
plot(merged$`Endothelial cell.x`, merged$`Endothelial cell.y`, xlab = "epic_tcga_firehouse_Endothelial_cells", ylab = "epic_tcga_pancancer_Endothelial_cells")
dev.off()

pdf("data-raw/benchmarking/tcga_firehouse_vs_pancancer_EPIC/Macrophage.pdf")
plot(merged$Macrophage.x, merged$Macrophage.y, xlab = "epic_tcga_firehouse_Macrophage", ylab = "epic_tcga_pancancer_Macrophage")
dev.off()

pdf("data-raw/benchmarking/tcga_firehouse_vs_pancancer_EPIC/NKcells.pdf")
plot(merged$`NK cell.x`, merged$`NK cell.y`, xlab = "epic_tcga_firehouse_NKcells", ylab = "epic_tcga_pancancer_NKcells")
dev.off()

pdf("data-raw/benchmarking/tcga_firehouse_vs_pancancer_EPIC/uncharacterized_cells.pdf")
plot(merged$`uncharacterized cell.x`, merged$`uncharacterized cell.y`, xlab = "epic_tcga_firehouse_uncharacterized_cells", ylab = "epic_tcga_pancancer_uncharacterized_cells")
dev.off()

# FORMATTING XCELL DATA

xcell_tcga_firehouse_t <- t(xcell_tcga_firehouse)
colnames(xcell_tcga_firehouse_t) <- xcell_tcga_firehouse_t[1, ]
xcell_tcga_firehouse_t <- xcell_tcga_firehouse_t[-1, ]
xcell_tcga_firehouse_t <- as.data.frame(xcell_tcga_firehouse_t)
xcell_tcga_firehouse_t[, 40] <- rownames(xcell_tcga_firehouse_t)

xcell_tcga_pancancer_t <- t(xcell_tcga_pancancer)
colnames(xcell_tcga_pancancer_t) <- xcell_tcga_pancancer_t[1, ]
xcell_tcga_pancancer_t <- xcell_tcga_pancancer_t[-1, ]
xcell_tcga_pancancer_t <- as.data.frame(xcell_tcga_pancancer_t)
xcell_tcga_pancancer_t[, 40] <- rownames(xcell_tcga_pancancer_t)

######################################################################################
# xcell - FIREHOUSE VS PANCANCER
######################################################################################

merged <- merge(xcell_tcga_firehouse_t, xcell_tcga_pancancer_t, by.x = "V40", by.y = "V40")

pdf("data-raw/benchmarking/tcga_firehouse_vs_pancancer_xcell/Bcells.pdf")
plot(merged$`B cell.x`, merged$`B cell.y`, xlab = "xcell_tcga_firehouse_Bcells", ylab = "xcell_tcga_pancancer_Bcells")
dev.off()
pdf("data-raw/benchmarking/tcga_firehouse_vs_pancancer_xcell/Myeloid_dendritic_cells.pdf")
plot(merged$`Myeloid dendritic cell activated.x`, merged$`Myeloid dendritic cell activated.y`, xlab = "xcell_tcga_firehouse_Myeloid_dendritic_cells_activated", ylab = "xcell_tcga_pancancer_Myeloid_dendritic_cells_activated")
dev.off()
pdf("data-raw/benchmarking/tcga_firehouse_vs_pancancer_xcell/Cancer_associated_fibroblasts.pdf")
plot(merged$`Cancer associated fibroblast.x`, merged$`Cancer associated fibroblast.y`, xlab = "xcell_tcga_firehouse_CAF", ylab = "xcell_tcga_pancancer_CAF")
dev.off()
pdf("data-raw/benchmarking/tcga_firehouse_vs_pancancer_xcell/T_cell_CD4+_memory.pdf")
plot(merged$`T cell CD4+ memory.x`, merged$`T cell CD4+ memory.y`, xlab = "xcell_tcga_firehouse_T_cell_CD4+_memory", ylab = "xcell_tcga_pancancer_T_cell_CD4+_memory")
dev.off()

pdf("data-raw/benchmarking/tcga_firehouse_vs_pancancer_xcell/T_cell_CD8+.pdf")
plot(merged$`T cell CD8+.x`, merged$`T cell CD8+.y`, xlab = "xcell_tcga_firehouse_T_cell_CD8+", ylab = "xcell_tcga_pancancer_T_cell_CD8+")
dev.off()
pdf("data-raw/benchmarking/tcga_firehouse_vs_pancancer_xcell/Endothelial_cells.pdf")
plot(merged$`Endothelial cell.x`, merged$`Endothelial cell.y`, xlab = "xcell_tcga_firehouse_Endothelial_cells", ylab = "xcell_tcga_pancancer_Endothelial_cells")
dev.off()
pdf("data-raw/benchmarking/tcga_firehouse_vs_pancancer_xcell/Macrophage.pdf")
plot(merged$Macrophage.x, merged$Macrophage.y, xlab = "xcell_tcga_firehouse_Macrophage", ylab = "xcell_tcga_pancancer_Macrophage")
dev.off()
pdf("data-raw/benchmarking/tcga_firehouse_vs_pancancer_xcell/NKcells.pdf")
plot(merged$`NK cell.x`, merged$`NK cell.y`, xlab = "xcell_tcga_firehouse_NKcells", ylab = "xcell_tcga_pancancer_NKcells")
dev.off()
pdf("data-raw/benchmarking/tcga_firehouse_vs_pancancer_xcell/myeloid_dendritic_cells.pdf")
plot(merged$`Myeloid dendritic cell.x`, merged$`Myeloid dendritic cell.y`, xlab = "xcell_tcga_firehouse_myeloid_dendritic_cells", ylab = "xcell_tcga_pancancer_myeloid_dendritic_cells")
dev.off()

pdf("data-raw/benchmarking/tcga_firehouse_vs_pancancer_xcell/Tcell_CD4+_naive.pdf")
plot(merged$`T cell CD4+ naive.x`, merged$`T cell CD4+ naive.y`, xlab = "xcell_tcga_firehouse_Tcell_CD4+_naive", ylab = "xcell_tcga_pancancer_Tcell_CD4+_naive")
dev.off()
pdf("data-raw/benchmarking/tcga_firehouse_vs_pancancer_xcell/Tcell_CD4+_non_regulatory.pdf")
plot(merged$`T cell CD4+ (non-regulatory).x`, merged$`T cell CD4+ (non-regulatory).y`, xlab = "xcell_tcga_firehouse_Tcell_CD4+_non_regulatory", ylab = "xcell_tcga_pancancer_Tcell_CD4+_non_regulatory")
dev.off()

# FORMATTING MCP DATA

mcp_tcga_firehouse_t <- t(mcp_tcga_firehouse)
colnames(mcp_tcga_firehouse_t) <- mcp_tcga_firehouse_t[1, ]
mcp_tcga_firehouse_t <- mcp_tcga_firehouse_t[-1, ]
mcp_tcga_firehouse_t <- as.data.frame(mcp_tcga_firehouse_t)
mcp_tcga_firehouse_t[, 12] <- rownames(mcp_tcga_firehouse_t)

mcp_tcga_pancancer_t <- t(mcp_tcga_pancancer)
colnames(mcp_tcga_pancancer_t) <- mcp_tcga_pancancer_t[1, ]
mcp_tcga_pancancer_t <- mcp_tcga_pancancer_t[-1, ]
mcp_tcga_pancancer_t <- as.data.frame(mcp_tcga_pancancer_t)
mcp_tcga_pancancer_t[, 12] <- rownames(mcp_tcga_pancancer_t)

######################################################################################
# MCP - FIREHOUSE VS PANCANCER
######################################################################################

merged <- merge(mcp_tcga_firehouse_t, mcp_tcga_pancancer_t, by.x = "V12", by.y = "V12")
pdf("data-raw/benchmarking/tcga_firehouse_vs_pancancer_mcp/Tcells.pdf")
plot(merged$`T cell.x`, merged$`T cell.y`, xlab = "mcp_tcga_firehouse_Tcells", ylab = "mcp_tcga_pancancer_Tcells")
dev.off()
pdf("data-raw/benchmarking/tcga_firehouse_vs_pancancer_mcp/Tcells_CD8+.pdf")
plot(merged$`T cell CD8+.x`, merged$`T cell CD8+.y`, xlab = "mcp_tcga_firehouse_Tcells_CD8+", ylab = "mcp_tcga_pancancer_Tcells_CD8+")
dev.off()
# pdf("data-raw/benchmarking/tcga_firehouse_vs_pancancer_mcp/T_cell_CD4+.pdf")
# plot(merged$`T cell CD4+.x`,merged$`T cell CD4+.y`,xlab="mcp_tcga_firehouse_T_cell_CD4+",ylab="mcp_tcga_pancancer_T_cell_CD4+")
# dev.off()
pdf("data-raw/benchmarking/tcga_firehouse_vs_pancancer_mcp/NKcells.pdf")
plot(merged$`NK cell.x`, merged$`NK cell.y`, xlab = "mcp_tcga_firehouse_NKcells", ylab = "mcp_tcga_pancancer_NKcells")
dev.off()
pdf("data-raw/benchmarking/tcga_firehouse_vs_pancancer_mcp/Bcells.pdf")
plot(merged$`B cell.x`, merged$`B cell.y`, xlab = "mcp_tcga_firehouse_Bcells", ylab = "mcp_tcga_pancancer_Bcells")
dev.off()
pdf("data-raw/benchmarking/tcga_firehouse_vs_pancancer_mcp/monocyte.pdf")
plot(merged$Monocyte.x, merged$Monocyte.y, xlab = "mcp_tcga_firehouse_monocyte", ylab = "mcp_tcga_pancancer_monocyte")
dev.off()
pdf("data-raw/benchmarking/tcga_firehouse_vs_pancancer_mcp/myeloid_dendritic_cells.pdf")
plot(merged$`Myeloid dendritic cell.x`, merged$`Myeloid dendritic cell.y`, xlab = "mcp_tcga_firehouse_myeloid_dendritic_cells", ylab = "mcp_tcga_pancancer_myeloid_dendritic_cells")
dev.off()
pdf("data-raw/benchmarking/tcga_firehouse_vs_pancancer_mcp/Endothelial_cells.pdf")
plot(merged$`Endothelial cell.x`, merged$`Endothelial cell.y`, xlab = "mcp_tcga_firehouse_Endothelial_cells", ylab = "mcp_tcga_pancancer_Endothelial_cells")
dev.off()

pdf("data-raw/benchmarking/tcga_firehouse_vs_pancancer_mcp/Neutrophils.pdf")
plot(merged$Neutrophil.x, merged$Neutrophil.y, xlab = "mcp_tcga_firehouse_neutrophils", ylab = "mcp_tcga_pancancer_neutrophils")
dev.off()


pdf("data-raw/benchmarking/tcga_firehouse_vs_pancancer_mcp/Macrophage.pdf")
plot(merged$Macrophage.x, merged$Macrophage.y, xlab = "mcp_tcga_firehouse_Macrophage", ylab = "mcp_tcga_pancancer_Macrophage")
dev.off()


pdf("data-raw/benchmarking/tcga_firehouse_vs_pancancer_mcp/Cancer_associated_fibroblasts.pdf")
plot(merged$`Cancer associated fibroblast.x`, merged$`Cancer associated fibroblast.y`, xlab = "mcp_tcga_firehouse_CAF", ylab = "mcp_tcga_pancancer_CAF")
dev.off()
