load("data-raw/gex_tcga.RData")
load("data-raw/osfgex_tcga.RData")

load("data-raw/gex_taylor.RData")
load("data-raw/gex_sun.RData")

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




