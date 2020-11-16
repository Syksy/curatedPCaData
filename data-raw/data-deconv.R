gex_tcga <- curatedPCaData::mae_tcga[["gex"]]

tcga_package_quantiseq <- immunedeconv::deconvolute(gex_tcga,"quantiseq")
tcga_package_epic <- immunedeconv::deconvolute(gex_tcga,"epic")
tcga_package_xcell <- immunedeconv::deconvolute(gex_tcga,"xcell")
tcga_package_mcp <- immunedeconv::deconvolute(gex_tcga,"mcp_counter")

write.csv(tcga_package_quantiseq, "data-raw/tcga_package_quantiseq.csv")
write.csv(tcga_package_epic, "data-raw/tcga_package_epic.csv")
write.csv(tcga_package_xcell, "data-raw/tcga_package_xcell.csv")
write.csv(tcga_package_mcp, "data-raw/tcga_package_mcp.csv")

gex_taylor <- curatedPCaData::mae_taylor[["gex"]]

taylor_package_quantiseq <- immunedeconv::deconvolute(gex_taylor,"quantiseq")
taylor_package_epic <- immunedeconv::deconvolute(gex_taylor,"epic")
taylor_package_xcell <- immunedeconv::deconvolute(gex_taylor,"xcell")
taylor_package_mcp <- immunedeconv::deconvolute(gex_taylor,"mcp_counter")

write.csv(taylor_package_quantiseq, "data-raw/taylor_package_quantiseq.csv")
write.csv(taylor_package_epic, "data-raw/taylor_package_epic.csv")
write.csv(taylor_package_xcell, "data-raw/taylor_package_xcell.csv")
write.csv(taylor_package_mcp, "data-raw/taylor_package_mcp.csv")


gex_sun <- curatedPCaData::mae_sun[["gex"]]

sun_package_quantiseq <- immunedeconv::deconvolute(gex_sun,"quantiseq")
sun_package_epic <- immunedeconv::deconvolute(gex_sun,"epic")
sun_package_xcell <- immunedeconv::deconvolute(gex_sun,"xcell")
sun_package_mcp <- immunedeconv::deconvolute(gex_sun,"mcp_counter")

write.csv(sun_package_quantiseq, "data-raw/sun_package_quantiseq.csv")
write.csv(sun_package_epic, "data-raw/sun_package_epic.csv")
write.csv(sun_package_xcell, "data-raw/sun_package_xcell.csv")
write.csv(sun_package_mcp, "data-raw/sun_package_mcp.csv")

