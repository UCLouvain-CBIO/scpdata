## Ensure resources are present

## code taken from 
## https://github.com/waldronlab/curatedTCGAData/blob/master/tests/testthat/test-resources.R

test_that("metadata numbers match ExperimentHub", {
    metaf <- system.file("extdata", "metadata.csv", 
                         package = "scpdata", mustWork = TRUE)
    meta <- read.csv(metaf, stringsAsFactors = FALSE)
    EHub <- query(ExperimentHub(), "scpdata")
    expect_identical(nrow(meta), length(EHub))
    expect_identical(meta$Title, EHub$title)
})

test_that("scpdata", {
    metaf <- system.file("extdata", "metadata.csv", 
                         package = "scpdata", mustWork = TRUE)
    meta <- read.csv(metaf, stringsAsFactors = FALSE)
    res <- scpdata()
    expect_identical(meta$Title, res$title)
})

test_that("all datasets are available", {
    res <- scpdata()
    for (dataset in res$title) {
        ds <- eval(call(dataset))
        expect_true(is(ds, "QFeatures"))
    }
})

