###
#
# Native 'R CMD check' tests run on the 'curatedPCaData'-package
# Any exceptions will count as a failure for 'R CMD check' run (notably, does not require 'RUnit' or 'testthat' packages for testing)
#
###

##
# Testing of getPCa main functionality
##

# Test retrieval of TCGA with all assays
# Get default fetching of a MAE object based on short id
methods::is(curatedPCaData::getPCa("tcga"), "MultiAssayExperiment")

# Test retrieval of Taylor with a pre-specified subset of assays
# Get fetching of an assay subset
methods::is(curatedPCaData::getPCa("taylor", assays = c("gex.rma", "cibersort", "scores")), "MultiAssayExperiment")

# Test a data fetch that should result in an error
# Test that an error is produced correctly for a study that does not exist
methods::is(try({curatedPCaData::getPCa("study_that_does_not_exist_or_is_misspelled", assays = c("foo", "bar"))}, silent=TRUE), "try-error")

# Test fetching of an assay that does not exist
# Test that an error is produced correctly for assays that do not exist
methods::is(try({curatedPCaData::getPCa("tcga", assays = "typo")}, silent=TRUE), "try-error")

##
# Testing of supporting summary functions etc
##

# Test fetching of study short ids and that the 19 studies originally available in Laajala et al. 2013 are retrieved correctly
# Tested function: curatedPCaData::getPCaStudies
studies <- curatedPCaData::getPCaStudies()
all(c("abida", "baca", "barbieri", "barwick", "chandran", "friedrich", "hieronymus", "icgcca", "igc", "kim", "kunderfranco", "ren", "sun", "taylor", "tcga", "true", "wallace", "wang", "weiner") %in% studies)

# Fetch MAE objects for further use
maes <- lapply(studies, FUN=\(id) { curatedPCaData::getPCa(id) })
names(maes) <- studies

# getPCaSummaryTable should summarize into a character matrix key instances and percentages for certain values for a given colData metadata variable
# Tested function: curatedPCaData::getPCaSummaryTable
inherits(curatedPCaData::getPCaSummaryTable(maes, var.name = "grade_group", vals=c("<=6", "3+4", "4+3", "7", ">=8")), "matrix")

# getPCaSummaryTable should summarize into a character matrix event counts and follow-up times for a Surv-like data
# Tested function: curatedPCaData::getPCaSummarySurv
inherits(curatedPCaData::getPCaSummarySurv(maes, event.name = "disease_specific_recurrence_status", time.name = "days_to_disease_specific_recurrence"), "matrix")


# getPCaSummarySamples should return a list of length 2; first element containing unique assay names and N counts in each study, and second element a matrix with GEX/CNA/MUT combinations for overlap
# Tested function: curatedPCaData::getPCaSummarySamples
inherits(curatedPCaData::getPCaSummarySamples(maes), "list")

# getPCaSummaryStudies should create a verbose character matrix depicting key characteristics for each study, such as sample counts, platforms, and special notes to be aware of
# Tested function: curatedPCaData::getPCaSummaryStudies, curatedPCaData::getPCaStudies
inherits(curatedPCaData::getPCaSummaryStudies(maes), "matrix")

