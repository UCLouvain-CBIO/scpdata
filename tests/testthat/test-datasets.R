## Ensure resources are present

## code taken from 
## https://github.com/waldronlab/curatedTCGAData/blob/master/tests/testthat/test-resources.R

test_that("metadata numbers match ExperimentHub", {
    assays_file <- system.file("extdata", 
                               "metadata.csv",
                               package = "scpdata",
                               mustWork = TRUE)
    metadataFile <- read.csv(assays_file, 
                             stringsAsFactors = FALSE)
    EHub <- query(ExperimentHub(), "scpdata")
    expect_equal(nrow(metadataFile), length(EHub))
    expect_equal(metadataFile$Title, EHub$title)
})