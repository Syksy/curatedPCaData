context("Testing that MAE object downloading from ExperimentHub/cache works and results in a functioning MAE object")

test_that("Testing getPCa and MAE construction", {
    expect_s4_class(curatedPCaData::getPCa("tcga"), "MultiAssayExperiment")
})

test_that("Testing getPCa and MAE construction", {
    expect_s4_class(curatedPCaData::getPCa("taylor"), "MultiAssayExperiment")
})
