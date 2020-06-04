#' TCGA gene expression. 
#'
#' A matrix containing the gene expression values for the 333 samples included 
#' in the 2015 Cell paper \url{https://www.ncbi.nlm.nih.gov/pubmed/26544944}
#'
#' @format A matrix with 333 rows and 36630 genes:
#' \describe{
#'   \item{genes}{raw gene expression values}
#'   ...
#' }
#' 
"tcga_gex"

#' TCGA copy number. 
#'
#' A matrix containing the copy number values for the 333 samples included 
#' in the 2015 Cell paper \url{https://www.ncbi.nlm.nih.gov/pubmed/26544944}
#'
#' @format A matrix with 333 rows and 36630 genes:
#' \describe{
#'   \item{genes}{copy number values}
#'   ...
#' }
#' 
"tcga_cna"

#' TCGA clinical/phenotype data. 
#'
#' A dataframe containing the clinical values for the 499 participant included 
#' in the provisional data set 
#'
#' @format A dataframe with 499 rows and 60 variable:
#' \describe{
#'   \item{study_name}{the study name that will link this information to the study meta-data}
#'   \item{patient_id}{a unique identifier for the patient}
#'   \item{sample_name}{primary sample identifier}
#'   \item{alt_sample_name}{if another identifier is used, for example in supplemental tables or GEO accession ids}
#'   \item{overall_survival_status}{binarized status of the patient, where 1 represents death and 0 represents no reported death}
#'   \item{days_to_overall_survival}{time to death or last follow-up in days}
#'   \item{age_at_initial_diagnosis}{in years}
#'   \item{year_diagnosis}{the year at which the patient was diagnoses with disease}
#'   \item{gleason_grade}{Gleason grade  (integer): total of two grades added together}
#'   \item{gleason_major}{Gleason grade  (integer) for the major site score}
#'   \item{gleason_minor}{Gleason grade  (integer) for the minor site score}
#'   \item{source_of_gleason}{source of where the pathogist performed Gleason grading}
#'   \item{grade_group}{separation of the gleason grades into groups that show different prognosis}
#'   \item{T_pathological}{pathological T stage (assessment made after surgery), based on tumor only}
#'   \item{T_substage_pathological}{pathological T substage (assessment made after surgery), based on tumor only}
#'   \item{T_clinical}{clinical T stage (at time of diagnosis) based on tumor only}
#'   \item{T_substage_clinical}{clinical T stage (at time of diagnosis) based on tumor only}
#'   \item{ERG_fusion_CNA}{presence of TMPRSS2:ERG gene fusion in prostate tumor as determined by copy number alteration analysis}
#'   \item{ERG_fusion_IHC}{presence of TMPRSS2:ERG gene fusion in prostate tumor as determined by immunohistochemistry}
#'   \item{ERG_fusion_GEX}{presence of TMPRSS2:ERG gene fusion in prostate tumor as determined by gene expression}
#'   \item{disease_specific_recurrence_status}{binarized status of the patient, where 1 represents an recurrence and 0 represents no reported recurrence}
#'   \item{days_to_disease_specific_recurrence}{time to recurrence or last follow-up in days}
#'   \item{metastasis_occurrence_status}{binarized status of the patient, where 1 represents a metastatic occurrence and 0 represents no reported metastatic occurrence}
#'   \item{days_to_metastatic_occurrence}{time to metastatic occurrence or last follow-up in days}
#'   \item{psa}{prostate specific angigen level at diagnosis}
#'   \item{race}{race}
#'   \item{smoking_status}{smoker (past or current) or non-smoker}
#'   \item{extraprostatic_extension}{spread of prostate cancer out of the prostate gland}
#'   \item{perineural_invasion}{cancer spreading to the space surrounding a nerve}
#'   \item{seminal_vesicle_invasion}{cancer has spread to the seminal vesicles}
#'   \item{angiolymphatic_invasion}{cancer has spread to blood vessels and lymph vessels}
#'   \item{androgen_ablation}{medical treatment to suppress or block the production of male sex hormones}
#'   \item{capsule}{status of the prostate capsule}
#'   \item{M_stage}{At the time of surgery. X: cannot evaluate distant metastasis, 0: there is no distant metastasis, 1: there is distant metastasis}
#'   \item{M_substage}{M1a: the cancer has spread to lymph nodes beyond the regional ones, M1b: the cancer has spread to_bone, M1c: the cancer has spread to other sites (regardless of bone involvement)}
#'   \item{other_patient}{a string that captures any additional patient information, features separated by bar (e.g. feature 1|feature 2|feature 3)}
#'   \item{sample_type}{type of tissue isolated from the patient and used for further -omic profiling}
#'   \item{genomic_alterations}{captures the list of reported alterations in the sample, in the format, gene:event, separated by a bar (eg. TP53:mutation|ETV1:fusion|PTEN:deletion)}
#'   \item{tumor_margins_positive}{Histologically altered cells in any surgical margins}
#'   \item{tissue_source}{the source of the sample}
#'   \item{metatstatic_site}{site where the metastatic sample was taken from}
#'   \item{microdissected}{microdissected or not}
#'   \item{frozen_ffpe}{frozen or FFPE}
#'   \item{other_feature}{other descriptions of the sample}
#'   \item{batch}{A way that describes a batch if a batch effect can be modeled (can be numeric or categorigal)}
#'   \item{other_sample}{a string that captures any additional sample information, features separated by bar (e.g. feature 1|feature 2|feature 3)}
#'   \item{tumor_purity_demix}{estimate of the tumor purity in the sample using the Demix method}
#'   \item{tumor_purity_absolute}{estimate of the tumor purity in the sample using Absolute}
#'   \item{zone_of_origin}{zone of origin assessed through tissue pathology}
#'   \item{mutational_signatures}{estimate of mutational signatures using deconstructSigs}
#'   \item{neoantigen_load}{estimate of mutational load using NetMHCPan}
#'   \item{AR_activity}{AR 20-gene signature}
#'   \item{N_stage}{regional lymph node status at the time of surgery.  X: cannot be measured, 0: no cancer in nearby lymph nodes, 1,2,3: the number and location of lymph nodes that contain cancer}
#'   \item{N_substage}{1a: the cancer has spread to lymph nodes beyond the regional ones, 1b: the cancer has spread to_bone, M1c: the cancer has spread to other sites (regardless of bone involvement)}
#'   \item{therapy_radiation_initial}{was radiation given as a primary therapy?}
#'   \item{therapy_radiation_salvage}{was radiation given after relapse from surgery?}
#'   \item{therapy_surgery_initial}{was surgery given as a primary therapy?}
#'   \item{therapy_hormonal_initial}{was hormonal therapy given as a primary therapy?}
#'   \item{other_treatment}{any other treatments}
#' }
#' 
"tcga_clinical"

#' Sun et al gene expression. 
#'
#' A matrix containing the gene expression values  \url{https://www.ncbi.nlm.nih.gov/pubmed/19343730}
#'
#' @format A matrix with 12057 rows and 79 samples:
#' \describe{
#'   \item{samples}{raw gene expression values for that sample}
#'   ...
#' }
#' 
"sun_gex"

#' Sun et al. clinical/phenotype data. 
#'
#' A dataframe containing the clinical values for the 79 participants 
#'
#' @format A dataframe with 79 rows and 62 variable:
#' \describe{
#'   \item{study_name}{the study name that will link this information to the study meta-data}
#'   \item{patient_id}{a unique identifier for the patient}
#'   \item{sample_name}{primary sample identifier}
#'   \item{alt_sample_name}{if another identifier is used, for example in supplemental tables or GEO accession ids}
#'   \item{overall_survival_status}{binarized status of the patient, where 1 represents death and 0 represents no reported death}
#'   \item{days_to_overall_survival}{time to death or last follow-up in days}
#'   \item{age_at_initial_diagnosis}{in years}
#'   \item{year_diagnosis}{the year at which the patient was diagnoses with disease}
#'   \item{gleason_grade}{Gleason grade  (integer): total of two grades added together}
#'   \item{gleason_major}{Gleason grade  (integer) for the major site score}
#'   \item{gleason_minor}{Gleason grade  (integer) for the minor site score}
#'   \item{source_of_gleason}{source of where the pathogist performed Gleason grading}
#'   \item{grade_group}{separation of the gleason grades into groups that show different prognosis}
#'   \item{T_pathological}{pathological T stage (assessment made after surgery), based on tumor only}
#'   \item{T_substage_pathological}{pathological T substage (assessment made after surgery), based on tumor only}
#'   \item{T_clinical}{clinical T stage (at time of diagnosis) based on tumor only}
#'   \item{T_substage_clinical}{clinical T stage (at time of diagnosis) based on tumor only}
#'   \item{ERG_fusion_CNA}{presence of TMPRSS2:ERG gene fusion in prostate tumor as determined by copy number alteration analysis}
#'   \item{ERG_fusion_IHC}{presence of TMPRSS2:ERG gene fusion in prostate tumor as determined by immunohistochemistry}
#'   \item{ERG_fusion_GEX}{presence of TMPRSS2:ERG gene fusion in prostate tumor as determined by gene expression}
#'   \item{disease_specific_recurrence_status}{binarized status of the patient, where 1 represents an recurrence and 0 represents no reported recurrence}
#'   \item{days_to_disease_specific_recurrence}{time to recurrence or last follow-up in days}
#'   \item{metastasis_occurrence_status}{binarized status of the patient, where 1 represents a metastatic occurrence and 0 represents no reported metastatic occurrence}
#'   \item{days_to_metastatic_occurrence}{time to metastatic occurrence or last follow-up in days}
#'   \item{psa}{prostate specific angigen level at diagnosis}
#'   \item{race}{race}
#'   \item{smoking_status}{smoker (past or current) or non-smoker}
#'   \item{extraprostatic_extension}{spread of prostate cancer out of the prostate gland}
#'   \item{perineural_invasion}{cancer spreading to the space surrounding a nerve}
#'   \item{seminal_vesicle_invasion}{cancer has spread to the seminal vesicles}
#'   \item{angiolymphatic_invasion}{cancer has spread to blood vessels and lymph vessels}
#'   \item{androgen_ablation}{medical treatment to suppress or block the production of male sex hormones}
#'   \item{capsule}{status of the prostate capsule}
#'   \item{M_stage}{At the time of surgery. X: cannot evaluate distant metastasis, 0: there is no distant metastasis, 1: there is distant metastasis}
#'   \item{M_substage}{M1a: the cancer has spread to lymph nodes beyond the regional ones, M1b: the cancer has spread to_bone, M1c: the cancer has spread to other sites (regardless of bone involvement)}
#'   \item{other_patient}{a string that captures any additional patient information, features separated by bar (e.g. feature 1|feature 2|feature 3)}
#'   \item{sample_type}{type of tissue isolated from the patient and used for further -omic profiling}
#'   \item{genomic_alterations}{captures the list of reported alterations in the sample, in the format, gene:event, separated by a bar (eg. TP53:mutation|ETV1:fusion|PTEN:deletion)}
#'   \item{tumor_margins_positive}{Histologically altered cells in any surgical margins}
#'   \item{tissue_source}{the source of the sample}
#'   \item{metatstatic_site}{site where the metastatic sample was taken from}
#'   \item{microdissected}{microdissected or not}
#'   \item{frozen_ffpe}{frozen or FFPE}
#'   \item{other_feature}{other descriptions of the sample}
#'   \item{batch}{A way that describes a batch if a batch effect can be modeled (can be numeric or categorigal)}
#'   \item{other_sample}{a string that captures any additional sample information, features separated by bar (e.g. feature 1|feature 2|feature 3)}
#'   \item{tumor_purity_demix}{estimate of the tumor purity in the sample using the Demix method}
#'   \item{tumor_purity_absolute}{estimate of the tumor purity in the sample using Absolute}
#'   \item{zone_of_origin}{zone of origin assessed through tissue pathology}
#'   \item{mutational_signatures}{estimate of mutational signatures using deconstructSigs}
#'   \item{neoantigen_load}{estimate of mutational load using NetMHCPan}
#'   \item{AR_activity}{AR 20-gene signature}
#'   \item{N_stage}{regional lymph node status at the time of surgery.  X: cannot be measured, 0: no cancer in nearby lymph nodes, 1,2,3: the number and location of lymph nodes that contain cancer}
#'   \item{N_substage}{1a: the cancer has spread to lymph nodes beyond the regional ones, 1b: the cancer has spread to_bone, M1c: the cancer has spread to other sites (regardless of bone involvement)}
#'   \item{therapy_radiation_initial}{was radiation given as a primary therapy?}
#'   \item{therapy_radiation_salvage}{was radiation given after relapse from surgery?}
#'   \item{therapy_surgery_initial}{was surgery given as a primary therapy?}
#'   \item{therapy_hormonal_initial}{was hormonal therapy given as a primary therapy?}
#'   \item{other_treatment}{any other treatments}
#' }
#' 
"sun_clinical"

#' Taylor et al. clinical/phenotype data. 
#'
#' A dataframe containing the clinical values for the 240 samples 
#'
#' @format A dataframe with 240 rows and 63 variable:
#' \describe{
#'   \item{study_name}{the study name that will link this information to the study meta-data}
#'   \item{patient_id}{a unique identifier for the patient}
#'   \item{sample_name}{primary sample identifier}
#'   \item{alt_sample_name}{if another identifier is used, for example in supplemental tables or GEO accession ids}
#'   \item{overall_survival_status}{binarized status of the patient, where 1 represents death and 0 represents no reported death}
#'   \item{days_to_overall_survival}{time to death or last follow-up in days}
#'   \item{age_at_initial_diagnosis}{in years}
#'   \item{year_diagnosis}{the year at which the patient was diagnoses with disease}
#'   \item{gleason_grade}{Gleason grade  (integer): total of two grades added together}
#'   \item{gleason_major}{Gleason grade  (integer) for the major site score}
#'   \item{gleason_minor}{Gleason grade  (integer) for the minor site score}
#'   \item{source_of_gleason}{source of where the pathogist performed Gleason grading}
#'   \item{grade_group}{separation of the gleason grades into groups that show different prognosis}
#'   \item{T_pathological}{pathological T stage (assessment made after surgery), based on tumor only}
#'   \item{T_substage_pathological}{pathological T substage (assessment made after surgery), based on tumor only}
#'   \item{T_clinical}{clinical T stage (at time of diagnosis) based on tumor only}
#'   \item{T_substage_clinical}{clinical T stage (at time of diagnosis) based on tumor only}
#'   \item{ERG_fusion_CNA}{presence of TMPRSS2:ERG gene fusion in prostate tumor as determined by copy number alteration analysis}
#'   \item{ERG_fusion_IHC}{presence of TMPRSS2:ERG gene fusion in prostate tumor as determined by immunohistochemistry}
#'   \item{ERG_fusion_GEX}{presence of TMPRSS2:ERG gene fusion in prostate tumor as determined by gene expression}
#'   \item{disease_specific_recurrence_status}{binarized status of the patient, where 1 represents an recurrence and 0 represents no reported recurrence}
#'   \item{days_to_disease_specific_recurrence}{time to recurrence or last follow-up in days}
#'   \item{metastasis_occurrence_status}{binarized status of the patient, where 1 represents a metastatic occurrence and 0 represents no reported metastatic occurrence}
#'   \item{days_to_metastatic_occurrence}{time to metastatic occurrence or last follow-up in days}
#'   \item{psa}{prostate specific angigen level at diagnosis}
#'   \item{race}{race}
#'   \item{smoking_status}{smoker (past or current) or non-smoker}
#'   \item{extraprostatic_extension}{spread of prostate cancer out of the prostate gland}
#'   \item{perineural_invasion}{cancer spreading to the space surrounding a nerve}
#'   \item{seminal_vesicle_invasion}{cancer has spread to the seminal vesicles}
#'   \item{angiolymphatic_invasion}{cancer has spread to blood vessels and lymph vessels}
#'   \item{androgen_ablation}{medical treatment to suppress or block the production of male sex hormones}
#'   \item{capsule}{status of the prostate capsule}
#'   \item{M_stage}{At the time of surgery. X: cannot evaluate distant metastasis, 0: there is no distant metastasis, 1: there is distant metastasis}
#'   \item{M_substage}{M1a: the cancer has spread to lymph nodes beyond the regional ones, M1b: the cancer has spread to_bone, M1c: the cancer has spread to other sites (regardless of bone involvement)}
#'   \item{other_patient}{a string that captures any additional patient information, features separated by bar (e.g. feature 1|feature 2|feature 3)}
#'   \item{sample_type}{type of tissue isolated from the patient and used for further -omic profiling}
#'   \item{genomic_alterations}{captures the list of reported alterations in the sample, in the format, gene:event, separated by a bar (eg. TP53:mutation|ETV1:fusion|PTEN:deletion)}
#'   \item{tumor_margins_positive}{Histologically altered cells in any surgical margins}
#'   \item{tissue_source}{the source of the sample}
#'   \item{metatstatic_site}{site where the metastatic sample was taken from}
#'   \item{microdissected}{microdissected or not}
#'   \item{frozen_ffpe}{frozen or FFPE}
#'   \item{other_feature}{other descriptions of the sample}
#'   \item{batch}{A way that describes a batch if a batch effect can be modeled (can be numeric or categorigal)}
#'   \item{other_sample}{a string that captures any additional sample information, features separated by bar (e.g. feature 1|feature 2|feature 3)}
#'   \item{tumor_purity_demix}{estimate of the tumor purity in the sample using the Demix method}
#'   \item{tumor_purity_absolute}{estimate of the tumor purity in the sample using Absolute}
#'   \item{zone_of_origin}{zone of origin assessed through tissue pathology}
#'   \item{mutational_signatures}{estimate of mutational signatures using deconstructSigs}
#'   \item{neoantigen_load}{estimate of mutational load using NetMHCPan}
#'   \item{AR_activity}{AR 20-gene signature}
#'   \item{N_stage}{regional lymph node status at the time of surgery.  X: cannot be measured, 0: no cancer in nearby lymph nodes, 1,2,3: the number and location of lymph nodes that contain cancer}
#'   \item{N_substage}{1a: the cancer has spread to lymph nodes beyond the regional ones, 1b: the cancer has spread to_bone, M1c: the cancer has spread to other sites (regardless of bone involvement)}
#'   \item{therapy_radiation_initial}{was radiation given as a primary therapy?}
#'   \item{therapy_radiation_salvage}{was radiation given after relapse from surgery?}
#'   \item{therapy_surgery_initial}{was surgery given as a primary therapy?}
#'   \item{therapy_hormonal_initial}{was hormonal therapy given as a primary therapy?}
#'   \item{other_treatment}{any other treatments}
#' }
#' 
"taylor_clinical"
