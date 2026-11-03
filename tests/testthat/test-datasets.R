## Ensure resources are present

## code taken from
## https://github.com/waldronlab/curatedTCGAData/blob/master/tests/testthat/test-resources.R

test_that("metadata numbers match ExperimentHub", {
    metaf <- system.file("extdata", "metadata.csv",
                         package = "scpdata", mustWork = TRUE)
    meta <- read.csv(metaf, stringsAsFactors = FALSE)
    EHub <- query(ExperimentHub(), "scpdata") ## same as mcols(scpdata())
    ## [2026-11-03 Tue] Changing the two expectations below, as they fail by
    ## default whenever a new experiment is added, i.e. regularly in devel. This
    ## unnecessarily increases maintenance burdens.
    ##
    ## expect_identical(nrow(meta), length(EHub))
    ## expect_identical(meta$Title, EHub$title)
    ##
    ## The new tests check inclusion of the remote EHub data in the local
    ## metadata.
    expect_true(nrow(meta) >= length(EHub))
    expect_true(all(EHub$title %in% meta$Title))
})

test_that("all datasets are available", {
    res <- scpdata()
    for (dataset in res$title) {
        ds <- eval(call(dataset))
        expect_true(is(ds, "QFeatures") | inherits(ds, "SummarizedExperiment"))
    }
})

test_that("leduc2022 throws a warning", {
    ## This data has been superseeded by leduc2022_plexDIA and leduc2022_pSCoPE.
    test <- expect_warning(leduc2022())
    expect_identical(length(test), 138L)
})