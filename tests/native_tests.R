###
#
# Native 'R CMD check' tests run on the 'curatedPCaData'-package
# Any exceptions will count as a failure for 'R CMD check' run (notably, does not require 'RUnit' or 'testthat' packages for testing)
#
###

# Test retrieval of TCGA with all assays
methods::is(curatedPCaData::getPCa("tcga"), "MultiAssayExperiment")

# Test retrieval of Taylor with a pre-specified subset of assays
methods::is(curatedPCaData::getPCa("taylor", assays = c("gex.rma", "cibersort", "scores")), "MultiAssayExperiment")

# Test a data fetch that should result in an error
methods::is(try({curatedPCaData::getPCa("study_that_does_not_exist_or_is_misspelled", assays = c("foo", "bar"))}, silent=TRUE), "try-error")

