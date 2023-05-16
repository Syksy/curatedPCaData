#############################################################################################################
# NEW CIBERSORT
#############################################################################################################
setwd("/Users/varsha/Desktop/gex_cibersortx/changed_gex_latest")

process_gex_cibersort <- function(gex) {
    gex <- as.data.frame(gex)
    gex[, ncol(gex) + 1] <- rownames(gex)
    rownames(gex) <- NULL
    gex <- gex[, c(ncol(gex), 1:(ncol(gex) - 1))]
}

gex_barwick <- curatedPCaData::mae_barwick[["gex.logq"]]
processed_gex_barwick <- process_gex_cibersort(gex_barwick)
write.table(processed_gex_barwick, "gex_barwick.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

gex_kunderfranco <- curatedPCaData::mae_kunderfranco[["gex.logr"]]
processed_gex_kunderfranco <- process_gex_cibersort(gex_kunderfranco)
write.table(processed_gex_kunderfranco, "gex_kunderfranco.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

gex_true <- curatedPCaData::mae_true[["gex.logr"]]
processed_gex_true <- process_gex_cibersort(gex_true)
write.table(processed_gex_true, "gex_true.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

gex_barbieri <- curatedPCaData::mae_barbieri[["gex.relz"]]
processed_gex_barbieri <- process_gex_cibersort(gex_barbieri)
write.table(processed_gex_barbieri, "gex_barbieri.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

gex_ren <- curatedPCaData::mae_ren[["gex.relz"]]
processed_gex_ren <- process_gex_cibersort(gex_ren)
write.table(processed_gex_ren, "gex_ren.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

gex_igc <- curatedPCaData::mae_igc[["gex.rma"]]
processed_gex_igc <- process_gex_cibersort(gex_igc)
write.table(processed_gex_igc, "gex_igc.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

gex_kim <- curatedPCaData::mae_kim[["gex.rma"]]
processed_gex_kim <- process_gex_cibersort(gex_kim)
write.table(processed_gex_kim, "gex_kim.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

gex_taylor <- curatedPCaData::mae_taylor[["gex.rma"]]
processed_gex_taylor <- process_gex_cibersort(gex_taylor)
write.table(processed_gex_taylor, "gex_taylor.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

gex_wallace <- curatedPCaData::mae_wallace[["gex.rma"]]
processed_gex_wallace <- process_gex_cibersort(gex_wallace)
write.table(processed_gex_wallace, "gex_wallace.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

gex_sun <- curatedPCaData::mae_sun[["gex.rma"]]
processed_gex_sun <- process_gex_cibersort(gex_sun)
write.table(processed_gex_sun, "gex_sun.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

gex_wang <- curatedPCaData::mae_wang[["gex.rma"]]
processed_gex_wang <- process_gex_cibersort(gex_wang)
write.table(processed_gex_wang, "gex_wang.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

gex_friedrich <- curatedPCaData::mae_friedrich[["gex.logq"]]
processed_gex_friedrich <- process_gex_cibersort(gex_friedrich)
write.table(processed_gex_friedrich, "gex_friedrich.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

gex_abida <- curatedPCaData::mae_abida[["gex.relz"]]
processed_gex_abida <- process_gex_cibersort(gex_abida)
write.table(processed_gex_abida, "gex_abida.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

gex_chandran <- curatedPCaData::mae_chandran[["gex.rma"]]
processed_gex_chandran <- process_gex_cibersort(gex_chandran)
write.table(processed_gex_chandran, "gex_chandran.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

gex_weiner <- curatedPCaData::mae_weiner[["gex.rma"]]
processed_gex_weiner <- process_gex_cibersort(gex_weiner)
write.table(processed_gex_weiner, "gex_weiner.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

gex_icgcca <- curatedPCaData::mae_icgcca[["gex.rma"]]
processed_gex_icgcca <- process_gex_cibersort(gex_icgcca)
write.table(processed_gex_icgcca, "gex_icgcca.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

gex_tcga <- curatedPCaData::mae_tcga[["gex.rsem.log"]]
processed_gex_tcga <- process_gex_cibersort(gex_tcga)
write.table(processed_gex_tcga, "gex_tcga.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
