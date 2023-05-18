#' Template for curatedPCaData clinical fields
#'
#' This data frame contains the clinical fields which were
#* the aim for extraction when gathering metadata for each dataset.
#' Variable types, ranges, and description is reported.
#'
#' @return A data.frame of the PCa template
#'
#' @docType data
#' @examples
#' data(template_prad)
#' head(template_prad)
"template_prad"

#' Abida et al. MAE-object
#'
#' MultiAssayExperiment object containing gene expression (gex), copy
#' number alteration (cna), mutations (mut) and immune cell estimates for Abida et al.
#'
#' @format A MAE object spanning 444 castrate resistant prostate cancer samples.
#' \describe{
#'    \item{cna.gistic}{matrix with 20291 rows and 444 columns, from GISTIC discretized copy number alteration calls.}
#'    \item{gex.relz}{matrix with 18971 rows and 266 columns, from z-score normalized expression in relative to paired normals.}
#'    \item{mut}{RaggedExperiment with 63184 rows and 444 columns, mutation data from cbioportal}
#'    \item{cibersort}{matrix with 22 rows and 261 columns, of cibersort based deconvolution data}
#'    \item{xcell}{matrix with 39 rows and 261 columns, of xcell based deconvolution data}
#'    \item{epic}{matrix with 8 rows and 261 columns, of epic based deconvolution data}
#'    \item{quantiseq}{matrix with 11 rows and 261 columns, of quantiseq based deconvolution data}
#'    \item{estimate}{data.frame with 4 rows and 261 columns, of cell types based on ESTIMATE method}
#'    \item{scores}{matrix with 4 rows and 261 columns, of risk scores and AR scores}
#'    \item{mcp}{mmatrix with 11 rows and 261 columns, of mcp-counter based deconvolution data}
#'
#' }
#' @details the clinical data refers to a sample of 444 tumors collected in 429 patients. The tissue was collected primarily at metastatic sites rather than from the prostate.
#' @references (\href{https://pubmed.ncbi.nlm.nih.gov/31061129/}{PubMed})
#' @source  \url{https://www.cbioportal.org/study/summary?id=prad_su2c_2019}
#' @return A MultiAssayExperiment corresponding to the study and its available omics.
#'
#' @rdname curatedPCaDatasets_abida
#' @examples
#' mae_abida <- getPCa("abida")
#' @name curatedPCaDatasets_abida
#' @keywords datasets
NULL

#' Baca et al. MAE-object
#'
#' MultiAssayExperiment object containing copy number alteration (cna) and mutations for Baca et al.
#'
#' @format A MAE object spanning 56 prostate cancer samples.
#' \describe{
#'   \item{cna.gistic}{matrix with 20124 rows and 56 columns, from GISTIC discretized copy number alteration calls.}
#'   \item{mut}{RaggedExperiment with 2734 rows and 57 columns, mutation data from cbioportal}
#' }
#' @source \url{https://www.cbioportal.org/study/summary?id=prad_broad_2013}
#' @return A MultiAssayExperiment corresponding to the study and its available omics.
#'
#' @rdname curatedPCaDatasets_baca
#' @docType data
#' @examples
#' mae_baca <- getPCa("baca")
#' @name curatedPCaDatasets_baca
#' @keywords datasets
NULL

#' Barbieri et al. MAE-object
#'
#' MultiAssayExperiment object containing gene expression (gex), copy
#' number alteration (cna), mutations (mut) and immune cell estimates for Barbieri et al.
#'
#' @format A MAE object spanning prostate adenocarcinomas from Barbieri et. al
#' \describe{
#'   \item{cna.gistic}{matrix with 20844 rows and 109 columns, from GISTIC discretized copy number alteration calls}
#'   \item{gex.relz}{matrix with 17917 rows and 31 columns, from z-score normalized expression in relative to paired normals.}
#'   \item{mut}{RaggedExperiment with 5737 rows and 112 columns, mutation data from cbioportal}
#'   \item{cibersort}{matrix with 22 rows and 31 columns, of cibersort based deconvolution data}
#'   \item{xcell}{matrix with 39 rows and 31 columns, of xcell based deconvolution data}
#'   \item{epic}{matrix with 8 rows and 31 columns, of epic based deconvolution data}
#'   \item{quantiseq}{matrix with 11 rows and 31 columns, of quantiseq based deconvolution data}
#'   \item{estimate}{data.frame with 4 rows and 31 columns, of cell types based on ESTIMATE method}
#'   \item{scores}{matrix with 4 rows and 31 columns, of risk scores and AR scores}
#'   \item{mcp}{matrix with 11 rows and 31 columns, of mcp-counter based deconvolution data}
#' }
#' @details The data comprises of primary localised prostate tumors from two cohorts, the Weill Cornell Medical College (WCMC; New York, NY), and the Uropath (Perth, Australia), which commercially provides banked urological tissues. None of the samples comes from patients who had received prior treatment for prostate cancer.
#' @references (\href{https://pubmed.ncbi.nlm.nih.gov/22610119/}{PubMed})
#' @source \url{https://www.cbioportal.org/study/summary?id=prad_broad}
#' @return A MultiAssayExperiment corresponding to the study and its available omics.
#'
#' @rdname curatedPCaDatasets_barbieri
#' @examples
#' mae_barbieri <- getPCa("barbieri")
#' @name curatedPCaDatasets_barbieri
#' @keywords datasets
NULL

#' Barwick et al. MAE-object
#'
#' MultiAssayExperiment object containing gene expression (gex) and immune cell estimates for Barwick et al.
#'
#' @format A MAE object spanning 146 prostate cancer samples.
#' \describe{
#'   \item{gex.logq}{matrix with 482 rows and 146 columns, for the log-quantile normalized gene expression data}
#'   \item{cibersort}{matrix with 22 rows and 139 columns, of cibersort based deconvolution data}
#'   \item{quantiseq}{matrix with 11 rows and 139 columns, of quantiseq based deconvolution data}
#'   \item{estimate}{data.frame with 4 rows and 139 columns, of cell types based on ESTIMATE method}
#'   \item{scores}{matrix with 3 rows and 139 columns, of risk scores and AR scores}
#'   \item{mcp}{matrix with 4 rows and 139 columns, of mcp-counter based deconvolution data}
#' }
#' @details Barwick et al. uses an older customized DASL array; therefore its gene coverage is lower, and many downstream methods fail due to lack of gene overlap.
#' @references (\href{https://pubmed.ncbi.nlm.nih.gov/22610119/}{PubMed})
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18655}
#' @return A MultiAssayExperiment corresponding to the study and its available omics.
#'
#' @rdname curatedPCaDatasets_barwick
#' @examples
#' mae_barwick <- getPCa("barwick")
#' @name curatedPCaDatasets_barwick
#' @keywords datasets
NULL

#' Chandran et al., Yu et al. MAE-object
#'
#' MultiAssayExperiment object containing gene expression (gex) and immune cell estimates for Chandran et al., Yu et al.
#'
#' @format An MAE object spanning 171 men with prostate cancer
#' \describe{
#'   \item{gex.rma}{matrix with 9007 rows and 171 columns, for the gene expression data}
#'   \item{cibersort}{matrix with 22 rows and 171 columns, of cibersort based deconvolution data}
#'   \item{xcell}{matrix with 39 rows and 171 columns, of xcell based deconvolution data}
#'   \item{epic}{matrix with 8 rows and 171 columns, of epic based deconvolution data}
#'   \item{quantiseq}{matrix with 11 rows and 171 columns, of quantiseq based deconvolution data}
#'   \item{estimate}{data.frame with 4 rows and 171 columns, of cell types based on ESTIMATE method}
#'   \item{scores}{matrix with 4 rows and 171 columns, of risk scores and AR scores}
#'   \item{mcp}{matrix with 11 rows and 171 columns, of mcp-counter based deconvolution data}
#' }
#' @references (\href{https://pubmed.ncbi.nlm.nih.gov/17430594/}{PubMed}) (\href{https://pubmed.ncbi.nlm.nih.gov/15254046/}{PubMed})
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6919}
#' @return A MultiAssayExperiment corresponding to the study and its available omics.
#'
#' @rdname curatedPCaDatasets_chandran
#' @examples
#' mae_chandran <- getPCa("chandran")
#' @name curatedPCaDatasets_chandran
#' @keywords datasets
NULL

#' Friedrich et al. MAE-object
#'
#' MultiAssayExperiment object containing gene expression (gex) and immune cell estimates for Friedrich et al.
#'
#' @format An MAE object spanning 255 men with prostate cancer
#' \describe{
#'   \item{gex.logq}{matrix with 23097 rows and 255 columns, for the log-quantile normalized gene expression data}
#'   \item{cibersort}{matrix with 22 rows and 255 columns, of cibersort based deconvolution data}
#'   \item{xcell}{matrix with 39 rows and 255 columns, of xcell based deconvolution data}
#'   \item{epic}{matrix with 8 rows and 255 columns, of epic based deconvolution data}
#'   \item{quantiseq}{matrix with 11 rows and 255 columns, of quantiseq based deconvolution data}
#'   \item{estimate}{data.frame with 4 rows and 255 columns, of cell types based on ESTIMATE method}
#'   \item{scores}{matrix with 4 rows and 255 columns, of risk scores and AR scores}
#'   \item{mcp}{matrix with 11 rows and 255 columns, of mcp-counter based deconvolution data}
#' }
#' @details The data comprises of 255 samples, with 164 primary tumors samples, 52 adjacent normal samples, and 39 benign prostate hyperplasia samples.  This dataset includes in its totality the 164 samples analysed in Kreutz et al. 2020.
#' @references (\href{https://pubmed.ncbi.nlm.nih.gov/32365858/}{PubMed}) (\href{https://pubmed.ncbi.nlm.nih.gov/32631745/}{PubMed})
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134051}
#' @return A MultiAssayExperiment corresponding to the study and its available omics.
#'
#' @rdname curatedPCaDatasets_friedrich
#' @examples
#' mae_friedrich <- getPCa("friedrich")
#' @name curatedPCaDatasets_friedrich
#' @keywords datasets
NULL

#' Hieronymus et al. MAE-object
#'
#' MultiAssayExperiment object containing copy number alteration (cna) for Hieronymus et al.
#'
#' @format A MAE object spanning 104 tumor samples
#' \describe{
#'   \item{cna.gistic}{matrix with 18026 rows and 104 columns, from GISTIC discretized copy number alteration calls}
#' }
#' @details The data comprises of 104 samples, for which are available clinical data and Copy Number Alteration, but no gene expression data -- thus no deconvolution results are available.
#' @references (\href{https://pubmed.ncbi.nlm.nih.gov/25024180/}{PubMed})
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54691}
#' @return A MultiAssayExperiment corresponding to the study and its available omics.
#'
#' @rdname curatedPCaDatasets_hieronymus
#' @examples
#' mae_hieronymus <- getPCa("hieronymus")
#' @name curatedPCaDatasets_hieronymus
#' @keywords datasets
NULL

#' ICGC CA MAE-object
#'
#' MultiAssayExperiment object containing gene expression (gex) and immune cell estimates for the ICGC CA (Canadian) cohort.
#'
#' @format An MAE object spanning 213 men with prostate cancer
#' \describe{
#'   \item{gex.rma}{matrix with 17208 rows and 213 columns, of gene expression data}
#'   \item{cibersort}{matrix with 22 rows and 213 columns, of cibersort based deconvolution data}
#'   \item{xcell}{matrix with 39 rows and 213 columns, of xcell based deconvolution data}
#'   \item{epic}{matrix with 8 rows and 213 columns, of epic based deconvolution data}
#'   \item{quantiseq}{matrix with 11 rows and 213 columns, of quantiseq based deconvolution data}
#'   \item{estimate}{data.frame with 4 rows and 213 columns, of cell types based on ESTIMATE method}
#'   \item{scores}{matrix with 4 rows and 213 columns, of risk scores and AR scores}
#'   \item{mcp}{matrix with 11 rows and 213 columns, of mcp-counter based deconvolution data}
#' }
#' @details The data refers to samples from the ICGC Canadian Prostate Cancer Genome Network (ICGC-PRAD-CA)
#' @references (\href{https://pubmed.ncbi.nlm.nih.gov/28068672/}{PubMed})
#' @source \url{https://dcc.icgc.org/projects/PRAD-CA}
#' @return A MultiAssayExperiment corresponding to the study and its available omics.
#'
#' @rdname curatedPCaDatasets_icgcca
#' @examples
#' mae_icgcca <- getPCa("icgcca")
#' @name curatedPCaDatasets_icgcca
#' @keywords datasets
NULL

#' International Genomics Consortium (IGC) MAE-object
#'
#' MultiAssayExperiment object containing gene expression (gex) and immune cell estimates from IGC's Expression Project for Oncology (expO) with a subset for prostate cancer
#'
#' @format An MAE object spanning 83 men with prostate cancer
#' \describe{
#'   \item{gex.rma}{matrix with 12798 rows and 83 columns, of gene expression data}
#'   \item{cibersort}{matrix with 22 rows and 83 columns, of cibersort based deconvolution data}
#'   \item{xcell}{matrix with 39 rows and 83 columns, of xcell based deconvolution data}
#'   \item{epic}{matrix with 8 rows and 83 columns, of epic based deconvolution data}
#'   \item{quantiseq}{matrix with 11 rows and 83 columns, of quantiseq based deconvolution data}
#'   \item{estimate}{data.frame with 4 rows and 83 columns, of cell types based on ESTIMATE method}
#'   \item{scores}{matrix with 4 rows and 83 columns, of risk scores and AR scores}
#'   \item{mcp}{matrix with 11 rows and 83 columns, of mcp-counter based deconvolution data}
#' }
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse2109}
#' @return A MultiAssayExperiment corresponding to the study and its available omics.
#'
#' @rdname curatedPCaDatasets_igc
#' @examples
#' mae_igc <- getPCa("igc")
#' @name curatedPCaDatasets_igc
#' @keywords datasets
NULL

#' Kim et al. MAE-object
#'
#' MultiAssayExperiment object containing gene expression (gex) and immune cell estimates from Kim et al.
#'
#' @format An MAE object spanning 266 men with prostate cancer
#' \describe{
#'   \item{gex.rma}{matrix with 17638 rows and 266 columns, of gene expression data}
#'   \item{cibersort}{matrix with 22 rows and 266 columns, of cibersort based deconvolution data}
#'   \item{xcell}{matrix with 39 rows and 266 columns, of xcell based deconvolution data}
#'   \item{epic}{matrix with 8 rows and 266 columns, of epic based deconvolution data}
#'   \item{quantiseq}{matrix with 11 rows and 266 columns, of quantiseq based deconvolution data}
#'   \item{estimate}{data.frame with 4 rows and 266 columns, of cell types based on ESTIMATE method}
#'   \item{scores}{matrix with 2 rows and 266 columns, of risk scores and AR scores}
#'   \item{mcp}{matrix with 11 rows and 266 columns, of mcp-counter based deconvolution data}
#' }
#' @details The dataset consists of 266 NCCN very low/low or favorable-intermediate risk PCa patients
#' @references (\href{https://pubmed.ncbi.nlm.nih.gov/30542054/}{PubMed})
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119616}
#' @return A MultiAssayExperiment corresponding to the study and its available omics.
#'
#' @rdname curatedPCaDatasets_kim
#' @examples
#' mae_kim <- getPCa("kim")
#' @name curatedPCaDatasets_kim
#' @keywords datasets
NULL

#' Kunderfranco et al. MAE-object
#'
#' MultiAssayExperiment object containing gene expression (gex) and immune cell estimates from Kunderfranco et al.
#'
#' @format An MAE object spanning 67 samples of normal prostate samples and prostate cancer samples
#' \describe{
#'   \item{gex.logr}{matrix with 16546 rows and 67 columns, of gene expression data}
#'   \item{cibersort}{matrix with 22 rows and 67 columns, of cibersort based deconvolution data}
#'   \item{xcell}{matrix with 39 rows and 67 columns, of xcell based deconvolution data}
#'   \item{epic}{matrix with 8 rows and 67 columns, of epic based deconvolution data}
#'   \item{quantiseq}{matrix with 11 rows and 67 columns, of quantiseq based deconvolution data}
#'   \item{estimate}{data.frame with 4 rows and 67 columns, of cell types based on ESTIMATE method}
#'   \item{scores}{matrix with 2 rows and 67 columns, of risk scores and AR scores}
#'   \item{mcp}{matrix with 11 rows and 67 columns, of mcp-counter based deconvolution data}
#' }
#' @details The data contains 14 disease free benign prostate hyperplasia samples and 53 prostate cancer samples.
#' @references (\href{https://pubmed.ncbi.nlm.nih.gov/20479932/}{PubMed})
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE14206}
#' @return A MultiAssayExperiment corresponding to the study and its available omics.
#'
#' @rdname curatedPCaDatasets_kunderfranco
#' @examples
#' mae_kunderfranco <- getPCa("kunderfranco")
#' @name curatedPCaDatasets_kunderfranco
#' @keywords datasets
NULL

#' Ren et al. MAE-object
#'
#' MultiAssayExperiment object containing gene expression (gex), mutations and immune cell estimates from Ren et al.
#'
#' @format An MAE object spanning 65 men with prostate cancer
#' \describe{
#'   \item{gex.relz}{matrix with 21046 rows and 65 columns, of gene expression data}
#'   \item{mut}{RaggedExperiment with 50625 rows and 65 columns, of mutation data from cbioportal}
#'   \item{cibersort}{matrix with 22 rows and 65 columns, of cibersort based deconvolution data}
#'   \item{xcell}{matrix with 39 rows and 65 columns, of xcell based deconvolution data}
#'   \item{epic}{matrix with 8 rows and 65 columns, of epic based deconvolution data}
#'   \item{quantiseq}{matrix with 11 rows and 65 columns, of quantiseq based deconvolution data}
#'   \item{estimate}{data.frame with 4 rows and 65 columns, of cell types based on ESTIMATE method}
#'   \item{scores}{matrix with 4 rows and 65 columns, of risk scores and AR scores}
#'   \item{mcp}{matrix with 11 rows and 65 columns, of mcp-counter based deconvolution data}
#' }
#' @details The data contains 14 disease free benign prostate hyperplasia smaples and 53 prostate cancer samples.
#' @references (\href{https://pubmed.ncbi.nlm.nih.gov/28927585/}{PubMed})
#' @source \url{https://www.cbioportal.org/study/summary?id=prad_eururol_2017}
#' @return A MultiAssayExperiment corresponding to the study and its available omics.
#'
#' @rdname curatedPCaDatasets_ren
#' @examples
#' mae_ren <- getPCa("ren")
#' @name curatedPCaDatasets_ren
#' @keywords datasets
NULL

#' Sun et al. MAE-object
#'
#' MultiAssayExperiment object containing gene expression (gex) and immune cell estimates from Sun et al.
#'
#' @format An MAE object spanning 79 men with prostate cancer
#' \describe{
#'   \item{gex.rma}{matrix with 12798 rows and 79 columns, of gene expression data}
#'   \item{cibersort}{matrix with 22 rows and 79 columns, of cibersort based deconvolution data}
#'   \item{xcell}{matrix with 39 rows and 79 columns, of xcell based deconvolution data}
#'   \item{epic}{matrix with 8 rows and 79 columns, of epic based deconvolution data}
#'   \item{quantiseq}{matrix with 11 rows and 79 columns, of quantiseq based deconvolution data}
#'   \item{estimate}{data.frame with 4 rows and 79 columns, of cell types based on ESTIMATE method}
#'   \item{scores}{matrix with 4 rows and 79 columns, of risk scores and AR scores}
#'   \item{mcp}{matrix with 11 rows and 79 columns, of mcp-counter based deconvolution data}
#' }
#'
#' @references (\href{https://pubmed.ncbi.nlm.nih.gov/19343730/}{PubMed})
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE25136}
#' @return A MultiAssayExperiment corresponding to the study and its available omics.
#'
#' @rdname curatedPCaDatasets_sun
#' @examples
#' mae_sun <- getPCa("sun")
#' @name curatedPCaDatasets_sun
#' @keywords datasets
NULL

#' Taylor et al. MAE-object
#'
#' MultiAssayExperiment object containing gene expression (gex), copy
#' number alteration (cna), mutations and immune cell estimates from Taylor et al.
#'
#' @format An MAE object spanning 218 men with prostate cancer
#' \describe{
#'   \item{cna.gistic}{matrix with 17832 rows and 194 columns, of gistic values for copy number alteration data}
#'   \item{cna.logr}{matrix with 18062 rows and 218 columns, of log-ratios for copy number alteration data}
#'   \item{gex.rma}{matrix with 17410 rows and 179 columns, of gene expression data}
#'   \item{mut}{RaggedExperiment with 90 rows and 43 columns, of mutation data}
#'   \item{cibersort}{matrix with 22 rows and 179 columns, of cibersort based deconvolution data}
#'   \item{xcell}{matrix with 39 rows and 179 columns, of xcell based deconvolution data}
#'   \item{epic}{matrix with 8 rows and 179 columns, of epic based deconvolution data}
#'   \item{quantiseq}{matrix with 11 rows and 179 columns, of quantiseq based deconvolution data}
#'   \item{estimate}{data.frame with 4 rows and 179 columns, of cell types based on ESTIMATE method}
#'   \item{scores}{matrix with 4 rows and 179 columns, of risk scores and AR scores}
#'   \item{mcp}{matrix with 11 rows and 179 columns, of mcp-counter based deconvolution data}
#' }
#'
#' @details Note that there is lack of overlap between the omics provided for each sample.
#' @references (\href{https://pubmed.ncbi.nlm.nih.gov/20579941/}{PubMed})
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE21035}
#' @return A MultiAssayExperiment corresponding to the study and its available omics.
#'
#' @rdname curatedPCaDatasets_taylor
#' @examples
#' mae_taylor <- getPCa("taylor")
#' @name curatedPCaDatasets_taylor
#' @keywords datasets
NULL

#' TCGA MAE-object
#'
#' MultiAssayExperiment object containing gene expression (gex), copy
#' number alteration (cna), mutations and immune cell estimates from GDC TCGA.
#'
#' @format An MAE object spanning 369 men with prostate cancer
#' \describe{
#'   \item{cna.gistic}{matrix with 23151 rows and 492 columns, from GISTIC discretized copy number alteration calls}
#'   \item{gex.rsem.log}{matrix with 19658 rows and 461 columns, of gene expression data}
#'   \item{mut}{RaggedExperiment with 30897 rows and 495 columns, of mutation data}
#'   \item{cibersort}{matrix with 22 rows and 550 columns, of cibersort based deconvolution data}
#'   \item{xcell}{matrix with 39 rows and 461 columns, of xcell based deconvolution data}
#'   \item{epic}{matrix with 8 rows and 461 columns, of epic based deconvolution data}
#'   \item{quantiseq}{matrix with 11 rows and 461 columns, of quantiseq based deconvolution data}
#'   \item{mcp}{matrix with 11 rows and 461 columns, of mcp-counter based deconvolution data}
#'   \item{estimate}{data.frame with 4 rows and 461 columns, of cell types based on ESTIMATE method}
#'   \item{scores}{matrix with 4 rows and 416 columns, of risk scores and AR scores}
#' }
#'
#' @details TCGA data was obtained from the latest GDC's XenaBrowser release.
#' @references (\href{https://pubmed.ncbi.nlm.nih.gov/26544944/}{PubMed})
#' @source \url{https://xenabrowser.net/datapages/?cohort=GDC\%20TCGA\%20Prostate\%20Cancer\%20(PRAD)}
#' @return A MultiAssayExperiment corresponding to the study and its available omics.
#'
#' @rdname curatedPCaDatasets_tcga
#' @examples
#' mae_tcga <- getPCa("tcga")
#' @name curatedPCaDatasets_tcga
#' @keywords datasets
NULL

#' True et al. MAE-object
#'
#' MultiAssayExperiment object containing gene expression (gex) and immune cell estimates from True et al.
#' @format An MAE object spanning 32 men with prostate cancer
#' \describe{
#'   \item{gex.logr}{matrix with 3615 rows and 32 columns, of gene expression data}
#'   \item{cibersort}{matrix with 22 rows and 32 columns, of cibersort based deconvolution data}
#'   \item{epic}{matrix with 8 rows and 32 columns, of epic based deconvolution data}
#'   \item{quantiseq}{matrix with 11 rows and 32 columns, of quantiseq based deconvolution data}
#'   \item{mcp}{matrix with 7 rows and 32 columns, of mcp-counter based deconvolution data}
#'   \item{scores}{matrix with 2 rows and 32 columns, of risk scores and AR scores}
#'   \item{estimate}{data.frame with 4 rows and 32 columns, of cell types based on ESTIMATE method}
#' }
#' @details The expression data for this dataset was produced by an older two colour chip that provides the relative expression levels for between tumor and normal.  Due to the chip type, there was a lack of gene overlap with some of the downstream methods, thus some downstream methods are missing.
# @references (\href{https://pubmed.ncbi.nlm.nih.gov/16829574/}{PubMed})
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5132}
#' @return A MultiAssayExperiment corresponding to the study and its available omics.
#'
#' @rdname curatedPCaDatasets_true
#' @examples
#' mae_true <- getPCa("true")
#' @name curatedPCaDatasets_true
#' @keywords datasets
NULL

#' Wallace et al. MAE-object
#'
#' MultiAssayExperiment object containing gene expression (gex) and immune cell estimates from Wallace et al.
#' @format An MAE object spanning 83 men with prostate cancer
#' \describe{
#'   \item{gex.rma}{matrix with 12783 rows and 89 columns, of gene expression data}
#'   \item{cibersort}{matrix with 22 rows and 83 columns, of cibersort based deconvolution data}
#'   \item{xcell}{matrix with 39 rows and 89 columns, of xcell based deconvolution data}
#'   \item{epic}{matrix with 8 rows and 89 columns, of epic based deconvolution data}
#'   \item{quantiseq}{matrix with 11 rows and 89 columns, of quantiseq based deconvolution data}
#'   \item{estimate}{data.frame with 4 rows and 89 columns, of cell types based on ESTIMATE method}
#'   \item{scores}{matrix with 4 rows and 89 columns, of risk scores and AR scores}
#'   \item{mcp}{matrix with 11 rows and 89 columns, of mcp-counter based deconvolution data}
#' }
#'
#' @references (\href{https://pubmed.ncbi.nlm.nih.gov/18245496/}{PubMed})
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6956}
#' @return A MultiAssayExperiment corresponding to the study and its available omics.
#'
#' @rdname curatedPCaDatasets_wallace
#' @examples
#' mae_wallace <- getPCa("wallace")
#' @name curatedPCaDatasets_wallace
#' @keywords datasets
NULL

#' Wang et al. MAE-object
#'
#' MultiAssayExperiment object containing gene expression (gex) and immune cell estimates from Wang et al.
#' @format An MAE object spanning 148 men with prostate cancer
#' \describe{
#'   \item{gex.rma}{matrix with 12798 rows and 148 columns, of gene expression data}
#'   \item{cibersort}{matrix with 22 rows and 148 columns, of cibersort based deconvolution data}
#'   \item{xcell}{matrix with 39 rows and 148 columns, of xcell based deconvolution data}
#'   \item{epic}{matrix with 8 rows and 148 columns, of epic based deconvolution data}
#'   \item{quantiseq}{matrix with 11 rows and 148 columns, of quantiseq based deconvolution data}
#'   \item{estimate}{data.frame with 4 rows and 148 columns, of cell types based on ESTIMATE method}
#'   \item{scores}{matrix with 4 rows and 148 columns, of risk scores and AR scores}
#'   \item{mcp}{matrix with 11 rows and 148 columns, of mcp-counter based deconvolution data}
#' }
#' @details 148 prostate samples, with various amounts of tumor, stroma, BPH and atrophic gland, were used for this study.
#' @references (\href{https://pubmed.ncbi.nlm.nih.gov/20663908/}{PubMed})
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE8218}
#' @return A MultiAssayExperiment corresponding to the study and its available omics.
#'
#' @rdname curatedPCaDatasets_wang
#' @examples
#' mae_wang <- getPCa("wang")
#' @name curatedPCaDatasets_wang
#' @keywords datasets
NULL

#' Weiner et al. MAE-object
#'
#' A MultiAssayExperiment object containing gene expression (gex) and immune cell estimates for Weiner et al.
#' @format An MAE spanning 838 prostate cancer samples of two cohorts
#' \describe{
#'  \item{gex.rma}{matrix of 17410 rows and 838 columns of gene expression data}
#'  \item{cibersort}{matrix with 22 rows and 838 columns, of cibersort based deconvolution data}
#'  \item{xcell}{matrix with 39 rows and 838 columns, the xcell deconvolution of the expression data}
#'  \item{epic}{matrix with 8 rows and 838 columns, the epic deconvolution of the expression data}
#'  \item{quantiseq}{matrix with 11 rows and 838 columns, the quantiseq deconvolution of the expression data}
#'  \item{estimate}{data.frame with 4 rows and 838 columns, of cell types based on ESTIMATE method}
#'  \item{scores}{matrix with 4 rows and 838 columns, of risk scores and AR scores}
#'  \item{mcp}{matrix with 11 rows and 838 columns, the mcp-counter deconvolution of the expression data}
#' }
#'
#' @references (\href{https://pubmed.ncbi.nlm.nih.gov/33568675/}{PubMed})
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157548}
#' @return A MultiAssayExperiment corresponding to the study and its available omics.
#'
#' @rdname curatedPCaDatasets_weiner
#' @examples
#' mae_weiner <- getPCa("weiner")
#' @name curatedPCaDatasets_weiner
#' @keywords datasets
NULL
