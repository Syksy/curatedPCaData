load("data-raw/gex_tcga.RData")
load("data-raw/osfgex_tcga.RData")

load("data-raw/gex_taylor.RData")
load("data-raw/gex_sun.RData")
load("data-raw/gex_barbieri.RData")
load("data-raw/gex_ren.RData")


osf_num <- as.matrix(osf)
mode(osf_num) <- "numeric"

# TCGA
tcga_package_quantiseq <- immunedeconv::deconvolute(gex_tcga,"quantiseq")
tcga_package_epic <- immunedeconv::deconvolute(gex_tcga,"epic")
tcga_package_xcell <- immunedeconv::deconvolute(gex_tcga,"xcell")
tcga_package_mcp <- immunedeconv::deconvolute(gex_tcga,"mcp_counter")

osf_quantiseq <- immunedeconv::deconvolute(osf_num,"quantiseq")
osf_epic <- immunedeconv::deconvolute(osf_num,"epic")
osf_xcell <- immunedeconv::deconvolute(osf_num,"xcell")
osf_mcp <- immunedeconv::deconvolute(osf_num,"mcp_counter")

# files have to be converted to txt to upload in the cibersort website
gex_tcga <- as.data.frame(gex_tcga)
gex_taylor <- as.data.frame(gex_taylor)
gex_sun <- as.data.frame(gex_sun)

data.table::setDT(gex_tcga,keep.rownames = "V1")
write.table(gex_tcga,"data-raw/tcga_gex.txt", sep="\t", row.names = F, col.names = T)

data.table::setDT(gex_taylor,keep.rownames = "V1")
write.table(gex_taylor,"data-raw/taylor_gex.txt", sep="\t", row.names = F, col.names = T)

data.table::setDT(gex_sun,keep.rownames = "V1")
write.table(gex_sun,"data-raw/sun_gex.txt", sep="\t", row.names = F, col.names = T)

# cibersortx results downloaded in csv format from "https://cibersortx.stanford.edu/" on 01/13/2021 at 12 pm EST using the gex matrices converted to txt files as shown above
osf_cibersortx <- rio::import("data-raw/CIBERSORTx_osf_Results.csv")
tcga_cibersortx <- rio::import("data-raw/CIBERSORTx_tcga_Results.csv")
taylor_cibersortx <- rio::import("data-raw/CIBERSORTx_taylor_Results.csv")


osf_cibersortx$Mixture <- gsub(osf_cibersortx$Mixture, pattern = "-", replacement = ".")
osf_cibersortx$Mixture <- paste(osf_cibersortx$Mixture, '01', sep='.')
osf_cibersortx_t <- t(osf_cibersortx)
colnames(osf_cibersortx_t) <- osf_cibersortx_t[1,]
osf_cibersortx_t <- osf_cibersortx_t[-1,]
save(osf_cibersortx_t,file="data-raw/osf_cibersortx_tcga.RData")

tcga_cibersortx_t <- t(tcga_cibersortx)
colnames(tcga_cibersortx_t) <- tcga_cibersortx_t[1,]
tcga_cibersortx_t <- tcga_cibersortx_t[-1,]
save(tcga_cibersortx_t,file="data-raw/cibersortx_tcga.RData")

tcga_package_quantiseq <- data.frame(tcga_package_quantiseq)
rownames(tcga_package_quantiseq) <- tcga_package_quantiseq[,1]
tcga_package_quantiseq <- tcga_package_quantiseq[,-1]
tcga_package_quantiseq <- as.matrix(tcga_package_quantiseq)
save(tcga_package_quantiseq,file="data-raw/quantiseq_tcga.RData")

osf_quantiseq <- data.frame(osf_quantiseq)
rownames(osf_quantiseq) <- osf_quantiseq[,1]
osf_quantiseq <- osf_quantiseq[,-1]
osf_quantiseq <- as.matrix(osf_quantiseq)
save(osf_quantiseq,file="data-raw/osf_quantiseq_tcga.RData")

tcga_package_epic <- data.frame(tcga_package_epic)
rownames(tcga_package_epic) <- tcga_package_epic[,1]
tcga_package_epic <- tcga_package_epic[,-1]
tcga_package_epic <- as.matrix(tcga_package_epic)
save(tcga_package_epic,file="data-raw/epic_tcga.RData")

osf_epic <- data.frame(osf_epic)
rownames(osf_epic) <- osf_epic[,1]
osf_epic <- osf_epic[,-1]
osf_epic <- as.matrix(osf_epic)
save(osf_epic,file="data-raw/osf_epic_tcga.RData")

tcga_package_xcell <- data.frame(tcga_package_xcell)
rownames(tcga_package_xcell) <- tcga_package_xcell[,1]
tcga_package_xcell <- tcga_package_xcell[,-1]
tcga_package_xcell <- as.matrix(tcga_package_xcell)
save(tcga_package_xcell,file="data-raw/xcell_tcga.RData")

osf_xcell <- data.frame(osf_xcell)
rownames(osf_xcell) <- osf_xcell[,1]
osf_xcell <- osf_xcell[,-1]
osf_xcell <- as.matrix(osf_xcell)
save(osf_xcell,file="data-raw/osf_xcell_tcga.RData")

tcga_package_mcp <- data.frame(tcga_package_mcp)
rownames(tcga_package_mcp) <- tcga_package_mcp[,1]
tcga_package_mcp <- tcga_package_mcp[,-1]
tcga_package_mcp <- as.matrix(tcga_package_mcp)
save(tcga_package_mcp,file="data-raw/mcp_tcga.RData")

osf_mcp <- data.frame(osf_mcp)
rownames(osf_mcp) <- osf_mcp[,1]
osf_mcp <- osf_mcp[,-1]
osf_mcp <- as.matrix(osf_mcp)
save(osf_mcp,file="data-raw/osf_mcp_tcga.RData")

mae_tcga <- create_mae(study_name = "TCGA")
usethis::use_data(mae_tcga, overwrite = TRUE)

# Taylor

taylor_package_quantiseq <- immunedeconv::deconvolute(gex_taylor,"quantiseq")
taylor_package_epic <- immunedeconv::deconvolute(gex_taylor,"epic")
taylor_package_xcell <- immunedeconv::deconvolute(gex_taylor,"xcell")
taylor_package_mcp <- immunedeconv::deconvolute(gex_taylor,"mcp_counter")

taylor_cibersortx_t <- t(taylor_cibersortx)
colnames(taylor_cibersortx_t) <- taylor_cibersortx_t[1,]
taylor_cibersortx_t <- taylor_cibersortx_t[-1,]
save(taylor_cibersortx_t,file="data-raw/cibersortx_taylor.RData")

taylor_package_quantiseq <- data.frame(taylor_package_quantiseq)
rownames(taylor_package_quantiseq) <- taylor_package_quantiseq[,1]
taylor_package_quantiseq <- taylor_package_quantiseq[,-1]
taylor_package_quantiseq <- as.matrix(taylor_package_quantiseq)
save(taylor_package_quantiseq,file="data-raw/quantiseq_taylor.RData")

taylor_package_epic <- data.frame(taylor_package_epic)
rownames(taylor_package_epic) <- taylor_package_epic[,1]
taylor_package_epic <- taylor_package_epic[,-1]
taylor_package_epic <- as.matrix(taylor_package_epic)
save(taylor_package_epic,file="data-raw/epic_taylor.RData")

taylor_package_xcell <- data.frame(taylor_package_xcell)
rownames(taylor_package_xcell) <- taylor_package_xcell[,1]
taylor_package_xcell <- taylor_package_xcell[,-1]
taylor_package_xcell <- as.matrix(taylor_package_xcell)
save(taylor_package_xcell,file="data-raw/xcell_taylor.RData")

taylor_package_mcp <- data.frame(taylor_package_mcp)
rownames(taylor_package_mcp) <- taylor_package_mcp[,1]
taylor_package_mcp <- taylor_package_mcp[,-1]
taylor_package_mcp <- as.matrix(taylor_package_mcp)
save(taylor_package_mcp,file="data-raw/mcp_taylor.RData")

mae_taylor <- create_mae(study_name = "taylor")
usethis::use_data(mae_taylor, internal = FALSE, overwrite = TRUE)

# Sun

sun_package_quantiseq <- immunedeconv::deconvolute(gex_sun,"quantiseq")
sun_package_epic <- immunedeconv::deconvolute(gex_sun,"epic")
sun_package_xcell <- immunedeconv::deconvolute(gex_sun,"xcell")
sun_package_mcp <- immunedeconv::deconvolute(gex_sun,"mcp_counter")

sun_package_quantiseq <- data.frame(sun_package_quantiseq)
rownames(sun_package_quantiseq) <- sun_package_quantiseq[,1]
sun_package_quantiseq <- sun_package_quantiseq[,-1]
sun_package_quantiseq <- as.matrix(sun_package_quantiseq)
save(sun_package_quantiseq,file="data-raw/quantiseq_sun.RData")

sun_package_epic <- data.frame(sun_package_epic)
rownames(sun_package_epic) <- sun_package_epic[,1]
sun_package_epic <- sun_package_epic[,-1]
sun_package_epic <- as.matrix(sun_package_epic)
save(sun_package_epic,file="data-raw/epic_sun.RData")

sun_package_xcell <- data.frame(sun_package_xcell)
rownames(sun_package_xcell) <- sun_package_xcell[,1]
sun_package_xcell <- sun_package_xcell[,-1]
sun_package_xcell <- as.matrix(sun_package_xcell)
save(sun_package_xcell,file="data-raw/xcell_sun.RData")

sun_package_mcp <- data.frame(sun_package_mcp)
rownames(sun_package_mcp) <- sun_package_mcp[,1]
sun_package_mcp <- sun_package_mcp[,-1]
sun_package_mcp <- as.matrix(sun_package_mcp)
save(sun_package_mcp,file="data-raw/mcp_sun.RData")

mae_sun <- create_mae(study_name = "Sun")
usethis::use_data(mae_sun, overwrite = TRUE)

# Barbieri Broad/Cornell data

Barbieri_package_quantiseq <- immunedeconv::deconvolute(gex_Barbieri,"quantiseq")
Barbieri_package_epic <- immunedeconv::deconvolute(gex_Barbieri,"epic")
Barbieri_package_xcell <- immunedeconv::deconvolute(gex_Barbieri,"xcell")
Barbieri_package_mcp <- immunedeconv::deconvolute(gex_Barbieri,"mcp_counter")

Barbieri_package_quantiseq <- data.frame(Barbieri_package_quantiseq)
rownames(Barbieri_package_quantiseq) <- Barbieri_package_quantiseq[,1]
Barbieri_package_quantiseq <- Barbieri_package_quantiseq[,-1]
Barbieri_package_quantiseq <- as.matrix(Barbieri_package_quantiseq)
save(Barbieri_package_quantiseq,file="data-raw/quantiseq_barbieri.RData")

Barbieri_package_epic <- data.frame(Barbieri_package_epic)
rownames(Barbieri_package_epic) <- Barbieri_package_epic[,1]
Barbieri_package_epic <- Barbieri_package_epic[,-1]
Barbieri_package_epic <- as.matrix(Barbieri_package_epic)
save(Barbieri_package_epic,file="data-raw/epic_barbieri.RData")

Barbieri_package_xcell <- data.frame(Barbieri_package_xcell)
rownames(Barbieri_package_xcell) <- Barbieri_package_xcell[,1]
Barbieri_package_xcell <- Barbieri_package_xcell[,-1]
Barbieri_package_xcell <- as.matrix(Barbieri_package_xcell)
save(Barbieri_package_xcell,file="data-raw/xcell_barbieri.RData")

Barbieri_package_mcp <- data.frame(Barbieri_package_mcp)
rownames(Barbieri_package_mcp) <- Barbieri_package_mcp[,1]
Barbieri_package_mcp <- Barbieri_package_mcp[,-1]
Barbieri_package_mcp <- as.matrix(Barbieri_package_mcp)
save(Barbieri_package_mcp,file="data-raw/mcp_barbieri.RData")

mae_barbieri <- create_mae(study_name = "Barbieri")
usethis::use_data(mae_barbieri, overwrite = TRUE)

# Ren Eururol

Ren_package_quantiseq <- immunedeconv::deconvolute(gex_Ren,"quantiseq")
Ren_package_epic <- immunedeconv::deconvolute(gex_Ren,"epic")
Ren_package_xcell <- immunedeconv::deconvolute(gex_Ren,"xcell")
Ren_package_mcp <- immunedeconv::deconvolute(gex_Ren,"mcp_counter")

Ren_package_quantiseq <- data.frame(Ren_package_quantiseq)
rownames(Ren_package_quantiseq) <- Ren_package_quantiseq[,1]
Ren_package_quantiseq <- Ren_package_quantiseq[,-1]
Ren_package_quantiseq <- as.matrix(Ren_package_quantiseq)
save(Ren_package_quantiseq,file="data-raw/quantiseq_ren.RData")

Ren_package_epic <- data.frame(Ren_package_epic)
rownames(Ren_package_epic) <- Ren_package_epic[,1]
Ren_package_epic <- Ren_package_epic[,-1]
Ren_package_epic <- as.matrix(Ren_package_epic)
save(Ren_package_epic,file="data-raw/epic_ren.RData")

Ren_package_xcell <- data.frame(Ren_package_xcell)
rownames(Ren_package_xcell) <- Ren_package_xcell[,1]
Ren_package_xcell <- Ren_package_xcell[,-1]
Ren_package_xcell <- as.matrix(Ren_package_xcell)
save(Ren_package_xcell,file="data-raw/xcell_ren.RData")

Ren_package_mcp <- data.frame(Ren_package_mcp)
rownames(Ren_package_mcp) <- Ren_package_mcp[,1]
Ren_package_mcp <- Ren_package_mcp[,-1]
Ren_package_mcp <- as.matrix(Ren_package_mcp)
save(Ren_package_mcp,file="data-raw/mcp_ren.RData")

mae_ren <- create_mae(study_name = "Ren")
usethis::use_data(mae_ren, overwrite = TRUE)



