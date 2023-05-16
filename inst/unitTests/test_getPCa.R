test_getPCa <- function(){
    testthat::expect_s4_class(curatedPCaData::getPCa("tcga"), "MultiAssayExperiment")
}    
