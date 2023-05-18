###
#
# Legacy roxygenization R-script for the lazyLoaded datasets
#
###


#' Template for curatedPCaData clinical fields
#'
#' This data frame contains the clinical fields which were the aim for extraction when gathering metadata for each dataset.
#' Variable types, ranges, and description is reported.
#'
#' @examples
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
#' @references Abida, W., Cyrta, J., Heller, G., et al. (2019). Genomic correlates of clinical outcome in advanced prostate cancer. Proceedings of the National Academy of Sciences of the United States of America, 116, 11428 - 11436.
#' (\href{https://pubmed.ncbi.nlm.nih.gov/31061129/}{PubMed})
#' @source  \url{https://www.cbioportal.org/study/summary?id=prad_su2c_2019}
#'
#' @examples
#' mae_abida
"mae_abida"

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
#'
#' @examples
#' mae_baca
"mae_baca"

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
#' @references Barbieri, C. E., Baca, S. C., Lawrence, M. S., et al. (2012). Exome sequencing identifies recurrent SPOP, FOXA1 and MED12 mutations in prostate cancer. Nature genetics, 44(6), 685–689. https://doi.org/10.1038/ng.2279
#' (\href{https://pubmed.ncbi.nlm.nih.gov/22610119/}{PubMed})
#' @source \url{https://www.cbioportal.org/study/summary?id=prad_broad}
#'
#' @examples
#' mae_barbieri
"mae_barbieri"

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
#' @references Barwick, B. G., Abramovitz, M., Kodani, M., Moreno, C. S., Nam, R., Tang, W., Bouzyk, M., Seth, A., & Leyland-Jones, B. (2010). Prostate cancer genes associated with TMPRSS2-ERG gene fusion and prognostic of biochemical recurrence in multiple cohorts. British journal of cancer, 102(3), 570–576. https://doi.org/10.1038/sj.bjc.6605519
#' (\href{https://pubmed.ncbi.nlm.nih.gov/22610119/}{PubMed})
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18655}
#'
#' @examples
#' mae_barwick
"mae_barwick"

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
#' @details .
#' @references Chandran, U. R., Ma, C., Dhir, R., et al. (2007). Gene expression profiles of prostate cancer reveal involvement of multiple molecular pathways in the metastatic process. BMC cancer, 7, 64. https://doi.org/10.1186/1471-2407-7-64
#' (\href{https://pubmed.ncbi.nlm.nih.gov/17430594/}{PubMed})
#' Yu, Y. P., Landsittel, D., Jing, L., et al. (2004). Gene expression alterations in prostate cancer predicting tumor aggression and preceding development of malignancy. Journal of clinical oncology : official journal of the American Society of Clinical Oncology, 22(14), 2790–2799. https://doi.org/10.1200/JCO.2004.05.158
#' (\href{https://pubmed.ncbi.nlm.nih.gov/15254046/}{PubMed})
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6919}
#'
#' @examples
#' mae_chandran
"mae_chandran"

#' Friedrich et al. MAE-object
#'
#' MultiAssayExperiment object containing gene expression (gex) and immune cell estimates for Friedrich et al.
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
#' @references Friedrich, M., Wiedemann, K., Reiche, K., et al. (2020). The Role of lncRNAs TAPIR-1 and -2 as Diagnostic Markers and Potential Therapeutic Targets in Prostate Cancer. Cancers, 12(5), 1122. https://doi.org/10.3390/cancers12051122
#' (\href{https://pubmed.ncbi.nlm.nih.gov/32365858/}{PubMed})
#' Kreuz, M., Otto, D. J., Fuessel, et al. (2020). ProstaTrend-A Multivariable Prognostic RNA Expression Score for Aggressive Prostate Cancer. European urology, 78(3), 452–459. https://doi.org/10.1016/j.eururo.2020.06.001
#' (\href{https://pubmed.ncbi.nlm.nih.gov/32631745/}{PubMed})
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134051}
#'
#' @examples
#' mae_friedrich
"mae_friedrich"

#' Hieronymus et al. MAE-object
#'
#' MultiAssayExperiment object containing copy number alteration (cna) for Hieronymus et al.
#'
#' @format A MAE object spanning 104 tumor samples
#' \describe{
#'   \item{cna.gistic}{matrix with 18026 rows and 104 columns, from GISTIC discretized copy number alteration calls}
#' }
#' @details The data comprises of 104 samples, for which are available clinical data and Copy Number Alteration, but no gene expression data -- thus no deconvolution results are available.
#' @references Hieronymus, H., Schultz, N., Gopalan, A., et al. (2014). Copy number alteration burden predicts prostate cancer relapse. Proceedings of the National Academy of Sciences of the United States of America, 111(30), 11139–11144. https://doi.org/10.1073/pnas.1411446111
#' (\href{https://pubmed.ncbi.nlm.nih.gov/25024180/}{PubMed})
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54691}
#'
#' @examples
#' mae_hieronymus
"mae_hieronymus"

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
#' @references Fraser, M., Sabelnykova, V.Y., Yamaguchi, T.N., et al. (2017). Genomic hallmarks of localized, non-indolent prostate cancer. Nature, 541, 359-364.
#'  (\href{https://pubmed.ncbi.nlm.nih.gov/28068672/}{PubMed})
#' @source \url{https://dcc.icgc.org/projects/PRAD-CA}
#'
#' @examples
#' mae_icgcca
"mae_icgcca"

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
#'
#' @examples
#' mae_igc
"mae_igc"

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
#' @references Kim HL, Li P, Huang HC, Deheshi S et al. Validation of the Decipher Test for predicting adverse pathology in candidates for prostate cancer active surveillance. Prostate Cancer Prostatic Dis 2019 Sep;22(3):399-405.
#' (\href{https://pubmed.ncbi.nlm.nih.gov/30542054/}{PubMed})
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119616}
#'
#' @examples
#' mae_kim
"mae_kim"

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
#' @details The data contains 14 disease free benign prostate hyperplasia smaples and 53 prostate cancer samples.
#' @references Kunderfranco, P., Mello-Grand, M., Cangemi, R., et al.  (2010). ETS transcription factors control transcription of EZH2 and epigenetic silencing of the tumor suppressor gene Nkx3.1 in prostate cancer. PloS one, 5(5), e10547. https://doi.org/10.1371/journal.pone.0010547
#'  (\href{https://pubmed.ncbi.nlm.nih.gov/20479932/}{PubMed})
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE14206}
#'
#' @examples
#' mae_kunderfranco
"mae_kunderfranco"

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
#' @references Ren, S., Wei, G. H., Liu, D., et al.  (2018). Whole-genome and Transcriptome Sequencing of Prostate Cancer Identify New Genetic Alterations Driving Disease Progression. European urology, 73(3), 322–339. https://doi.org/10.1016/j.eururo.2017.08.027
#'  (\href{https://pubmed.ncbi.nlm.nih.gov/28927585/}{PubMed})
#' @source \url{https://www.cbioportal.org/study/summary?id=prad_eururol_2017}
#'
#' @examples
#' mae_ren
"mae_ren"

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
#' @references Sun, Y., & Goodison, S. (2009). Optimizing molecular signatures for predicting prostate cancer recurrence. The Prostate, 69(10), 1119–1127. https://doi.org/10.1002/pros.20961
#'  (\href{https://pubmed.ncbi.nlm.nih.gov/19343730/}{PubMed})
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE25136}
#'
#' @examples
#' mae_sun
"mae_sun"

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
#' @references Taylor, B. S., Schultz, N., Hieronymus, H., et al. (2010). Integrative genomic profiling of human prostate cancer. Cancer cell, 18(1), 11–22. https://doi.org/10.1016/j.ccr.2010.05.026
#'  (\href{https://pubmed.ncbi.nlm.nih.gov/20579941/}{PubMed})
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE21035}
#'
#' @examples
#' mae_taylor
"mae_taylor"

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
#' @references Cancer Genome Atlas Research Network (2015). The Molecular Taxonomy of Primary Prostate Cancer. Cell, 163(4), 1011-1025. https://doi.org/10.1016/j.cell.2015.10.025
#'  (\href{https://pubmed.ncbi.nlm.nih.gov/26544944/}{PubMed})
#' @source \url{https://xenabrowser.net/datapages/?cohort=GDC\%20TCGA\%20Prostate\%20Cancer\%20(PRAD)}
#'
#' @examples
#' mae_tcga
"mae_tcga"

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
#' @references True, L., Coleman, I., Hawley, S., et al. (2006). A molecular correlate to the Gleason grading system for prostate adenocarcinoma. Proceedings of the National Academy of Sciences of the United States of America, 103(29), 10991–10996. https://doi.org/10.1073/pnas.0603678103
#'  (\href{https://pubmed.ncbi.nlm.nih.gov/16829574/}{PubMed})
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5132}
#'
#' @examples
#' mae_true
"mae_true"

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
#' @references Wallace, T. A., Prueitt, R. L., Yi, M., Howe, T. M., Gillespie, J. W., Yfantis, H. G., Stephens, R. M., Caporaso, N. E., Loffredo, C. A., & Ambs, S. (2008). Tumor immunobiological differences in prostate cancer between African-American and European-American men. Cancer research, 68(3), 927–936. https://doi.org/10.1158/0008-5472.CAN-07-2608
#'  (\href{https://pubmed.ncbi.nlm.nih.gov/18245496/}{PubMed})
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6956}
#'
#' @examples
#' mae_wallace
"mae_wallace"

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
#' @references Wang Y, Xia XQ, Jia Z, Sawyers A et al. In silico estimates of tissue components in surgical samples based on expression profiling data. Cancer Res 2010 Aug 15;70(16):6448-55.
#'  (\href{https://pubmed.ncbi.nlm.nih.gov/20663908/}{PubMed})
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE8218}
#'
#' @examples
#' mae_wang
"mae_wang"


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
#' @references Weiner, A.B., Vidotto, T., Liu, Y. et al. Plasma cells are enriched in localized prostate cancer in Black men and are associated with improved outcomes. Nat Commun 12, 935 (2021). https://doi.org/10.1038/s41467-021-21245-w
#' (\href{https://pubmed.ncbi.nlm.nih.gov/33568675/}{PubMed})
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157548}
#'
#' @examples
#' mae_weiner
"mae_weiner"
