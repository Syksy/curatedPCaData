load("data-raw/gex_tcga.RData")
load("data-raw/gex_taylor.RData")
load("data-raw/gex_sun.RData")

# TCGA
tcga_package_quantiseq <- immunedeconv::deconvolute(gex_tcga,"quantiseq")
tcga_package_epic <- immunedeconv::deconvolute(gex_tcga,"epic")
tcga_package_xcell <- immunedeconv::deconvolute(gex_tcga,"xcell")
tcga_package_mcp <- immunedeconv::deconvolute(gex_tcga,"mcp_counter")

tcga_package_quantiseq_t <- t(tcga_package_quantiseq)
colnames(tcga_package_quantiseq_t) <- tcga_package_quantiseq_t[1,]
tcga_package_quantiseq_t <- tcga_package_quantiseq_t[-1,]
tcga_package_quantiseq_t <- data.frame(tcga_package_quantiseq_t)
colnames(tcga_package_quantiseq_t) <- paste("quantiseq", colnames(tcga_package_quantiseq_t), sep = "_")
data.table::setDT(tcga_package_quantiseq_t,keep.rownames = "sample_name")
#tcga_package_quantiseq_t$sample_name <-gsub("\\.01*","",tcga_package_quantiseq_t$sample_name)



tcga_package_epic_t <- t(tcga_package_epic)
colnames(tcga_package_epic_t) <- tcga_package_epic_t[1,]
tcga_package_epic_t <- tcga_package_epic_t[-1,]
tcga_package_epic_t <- data.frame(tcga_package_epic_t)
colnames(tcga_package_epic_t) <- paste("epic", colnames(tcga_package_epic_t), sep = "_")
data.table::setDT(tcga_package_epic_t,keep.rownames = "sample_name")
#tcga_package_epic_t$sample_name <-gsub("\\.01*","",tcga_package_epic_t$sample_name)


tcga_package_xcell_t <- t(tcga_package_xcell)
colnames(tcga_package_xcell_t) <- tcga_package_xcell_t[1,]
tcga_package_xcell_t <- tcga_package_xcell_t[-1,]
tcga_package_xcell_t <- data.frame(tcga_package_xcell_t)
colnames(tcga_package_xcell_t) <- paste("xcell", colnames(tcga_package_xcell_t), sep = "_")
data.table::setDT(tcga_package_xcell_t,keep.rownames = "sample_name")
#tcga_package_xcell_t$sample_name <-gsub("\\.01*","",tcga_package_xcell_t$sample_name)

tcga_package_mcp_t <- t(tcga_package_mcp)
colnames(tcga_package_mcp_t) <- tcga_package_mcp_t[1,]
tcga_package_mcp_t <- tcga_package_mcp_t[-1,]
tcga_package_mcp_t <- data.frame(tcga_package_mcp_t)
colnames(tcga_package_mcp_t) <- paste("mcp", colnames(tcga_package_mcp_t), sep = "_")
data.table::setDT(tcga_package_mcp_t,keep.rownames = "sample_name")
#tcga_package_mcp_t$sample_name <-gsub("\\.01*","",tcga_package_mcp_t$sample_name)

merged_quantiseq_epic_tcga <- merge(tcga_package_quantiseq_t,tcga_package_epic_t,by="sample_name")
merged_xcell_mcp_tcga <- merge(tcga_package_xcell_t,tcga_package_mcp_t,by="sample_name")
deconv_tcga <- merge(merged_quantiseq_epic_tcga,merged_xcell_mcp_tcga,by="sample_name")
deconv_tcga <- data.frame(deconv_tcga, row.names = 1)
deconv_tcga <- t(deconv_tcga)

#deconv_tcga1 <- as.matrix(deconv_tcga)


save(deconv_tcga,file="data-raw/deconv_tcga.RData")

# Taylor
taylor_package_quantiseq <- immunedeconv::deconvolute(gex_taylor,"quantiseq")
taylor_package_epic <- immunedeconv::deconvolute(gex_taylor,"epic")
taylor_package_xcell <- immunedeconv::deconvolute(gex_taylor,"xcell")
taylor_package_mcp <- immunedeconv::deconvolute(gex_taylor,"mcp_counter")

taylor_package_quantiseq_t <- t(taylor_package_quantiseq)
colnames(taylor_package_quantiseq_t) <- taylor_package_quantiseq_t[1,]
taylor_package_quantiseq_t <- taylor_package_quantiseq_t[-1,]
taylor_package_quantiseq_t <- data.frame(taylor_package_quantiseq_t)
colnames(taylor_package_quantiseq_t) <- paste("quantiseq", colnames(taylor_package_quantiseq_t), sep = "_")
data.table::setDT(taylor_package_quantiseq_t,keep.rownames = "sample_name")

taylor_package_epic_t <- t(taylor_package_epic)
colnames(taylor_package_epic_t) <- taylor_package_epic_t[1,]
taylor_package_epic_t <- taylor_package_epic_t[-1,]
taylor_package_epic_t <- data.frame(taylor_package_epic_t)
colnames(taylor_package_epic_t) <- paste("epic", colnames(taylor_package_epic_t), sep = "_")
data.table::setDT(taylor_package_epic_t,keep.rownames = "sample_name")

taylor_package_xcell_t <- t(taylor_package_xcell)
colnames(taylor_package_xcell_t) <- taylor_package_xcell_t[1,]
taylor_package_xcell_t <- taylor_package_xcell_t[-1,]
taylor_package_xcell_t <- data.frame(taylor_package_xcell_t)
colnames(taylor_package_xcell_t) <- paste("xcell", colnames(taylor_package_xcell_t), sep = "_")
data.table::setDT(taylor_package_xcell_t,keep.rownames = "sample_name")

taylor_package_mcp_t <- t(taylor_package_mcp)
colnames(taylor_package_mcp_t) <- taylor_package_mcp_t[1,]
taylor_package_mcp_t <- taylor_package_mcp_t[-1,]
taylor_package_mcp_t <- data.frame(taylor_package_mcp_t)
colnames(taylor_package_mcp_t) <- paste("mcp", colnames(taylor_package_mcp_t), sep = "_")
data.table::setDT(taylor_package_mcp_t,keep.rownames = "sample_name")

merged_quantiseq_epic_taylor <- merge(taylor_package_quantiseq_t,taylor_package_epic_t,by="sample_name")
merged_xcell_mcp_taylor <- merge(taylor_package_xcell_t,taylor_package_mcp_t,by="sample_name")
deconv_taylor <- merge(merged_quantiseq_epic_taylor,merged_xcell_mcp_taylor,by="sample_name")

deconv_taylor$sample_name <- gsub("\\_.*", "", deconv_taylor$sample_name)
gsub("\\..*","",a)

deconv_taylor <- data.frame(deconv_taylor, row.names = 1)

deconv_taylor <- t(deconv_taylor)

#deconv_tcga1 <- as.matrix(deconv_tcga)

save(deconv_taylor,file="data-raw/deconv_taylor.RData")

# Sun

sun_package_quantiseq <- immunedeconv::deconvolute(gex_sun,"quantiseq")
sun_package_epic <- immunedeconv::deconvolute(gex_sun,"epic")
sun_package_xcell <- immunedeconv::deconvolute(gex_sun,"xcell")
sun_package_mcp <- immunedeconv::deconvolute(gex_sun,"mcp_counter")

sun_package_quantiseq_t <- t(sun_package_quantiseq)
colnames(sun_package_quantiseq_t) <- sun_package_quantiseq_t[1,]
sun_package_quantiseq_t <- sun_package_quantiseq_t[-1,]
sun_package_quantiseq_t <- data.frame(sun_package_quantiseq_t)
colnames(sun_package_quantiseq_t) <- paste("quantiseq", colnames(sun_package_quantiseq_t), sep = "_")
data.table::setDT(sun_package_quantiseq_t,keep.rownames = "V1")

sun_package_epic_t <- t(sun_package_epic)
colnames(sun_package_epic_t) <- sun_package_epic_t[1,]
sun_package_epic_t <- sun_package_epic_t[-1,]
sun_package_epic_t <- data.frame(sun_package_epic_t)
colnames(sun_package_epic_t) <- paste("epic", colnames(sun_package_epic_t), sep = "_")
data.table::setDT(sun_package_epic_t,keep.rownames = "V1")

sun_package_xcell_t <- t(sun_package_xcell)
colnames(sun_package_xcell_t) <- sun_package_xcell_t[1,]
sun_package_xcell_t <- sun_package_xcell_t[-1,]
sun_package_xcell_t <- data.frame(sun_package_xcell_t)
colnames(sun_package_xcell_t) <- paste("xcell", colnames(sun_package_xcell_t), sep = "_")
data.table::setDT(sun_package_xcell_t,keep.rownames = "V1")

sun_package_mcp_t <- t(sun_package_mcp)
colnames(sun_package_mcp_t) <- sun_package_mcp_t[1,]
sun_package_mcp_t <- sun_package_mcp_t[-1,]
sun_package_mcp_t <- data.frame(sun_package_mcp_t)
colnames(sun_package_mcp_t) <- paste("mcp", colnames(sun_package_mcp_t), sep = "_")
data.table::setDT(sun_package_mcp_t,keep.rownames = "V1")

merged_quantiseq_epic_sun<- merge(sun_package_quantiseq_t,sun_package_epic_t,by="V1")
merged_xcell_mcp_sun <- merge(sun_package_xcell_t,sun_package_mcp_t,by="V1")
deconv_taylor <- merge(merged_quantiseq_epic_sun,merged_xcell_mcp_sun,by="V1")

save(deconv_sun,file="data-raw/deconv_sun.RData")




