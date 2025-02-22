stopifnot(packageVersion("MsDataHub") >= "1.7.1")
stopifnot(packageVersion("QFeatures") >= "1.17.2")

library(tidyverse)
library(QFeatures)

ai2025aTab <- read_tsv(MsDataHub::Ai2025_aCMs_report.tsv())


tab <- tibble(File.Name = unique(ai2025aTab[[1]])) |>
    mutate(Sample = sub("^.+CM_PROJECT\\\\", "", File.Name)) |>
    mutate(Sample = sub("\\\\", "_", Sample)) |>
    mutate(Date = ymd(as.integer(substring(Sample, 1, 6)))) |>
    mutate(Subject = sub("^.+_(Subject[0-9])_.+$", "\\1", Sample)) |>
    mutate(PlateWell = sub("^.+_([A-Z][0-9]+)_.+$", "\\1", Sample)) |>
    mutate(PlateRow = gsub("[0-9]*", "", PlateWell)) |>
    mutate(PlateColumn = gsub("[A-Z]*", "", PlateWell)) |>
    mutate(PlateColumn = factor(PlateColumn, levels = 1:24)) |>
    mutate(Position = sub("^.+_([A-Z]+[0-9]+)_1_[0-9]+\\.d$", "\\1", Sample)) |>
    mutate(FileIndex = sub("^.+_1_([0-9]+)\\.d$", "\\1", Sample)) |>
    mutate(FileIndex = factor(FileIndex, levels = sort(unique(FileIndex))))

## Extracting the heart locations
tab$HeartLocation <- NA
ExpectedLocations <- c("Lvendo", "Lvepi", "Lvmid", "RV", "sytox")
for (i in 1:5) {
  loc <- ExpectedLocations[i]
  tab$HeartLocation[grep(loc, tab$File.Name, ignore.case = TRUE)] <- loc
}
tab$runCol <- tab$File.Name

ai2025a <- readSCPfromDIANN(ai2025aTab,
                         colData = DataFrame(tab),
                         fnames = "Precursor.Id")

## Setting the names of the QFeatures objects
ai2025a <- MultiAssayExperiment::renamePrimary(ai2025a, ai2025a$Sample)
names(ai2025a) <- ai2025a$Sample


#########################################################
## Processing from the vignette
acms <- zeroIsNA(ai2025a, names(ai2025a))
acms <- filterFeatures(acms, ~ abs(RT - Predicted.RT) < 0.18) |>
    filterFeatures(~ !grepl(";", Protein.Names) &
                       Proteotypic == 1)


acms$MedianIntensity <- sapply(names(acms), function(i) {
    median(log2(assay(acms[[i]])), na.rm = TRUE)
})
acms$TotalIds <- nrows(acms)

acms <- joinAssays(acms, i = names(acms), name = "precursors")

acms <- logTransform(acms, "precursors", "precursors_log")

acms <- aggregateFeatures(acms,
                          i = "precursors_log",
                          name = "peptides",
                          fcol = "Modified.Sequence",
                          fun = colMedians,
                          na.rm = TRUE)

acms <- aggregateFeatures(acms,
                          i = "peptides",
                          name = "proteins",
                          fcol = "Protein.Ids",
                          fun = colMedians,
                          na.rm = TRUE)

#########################################################
## Statistical modelling

sce <- getWithColData(acms, "precursors_log")
sce <- sce[, sce$HeartLocation != "sytox"]

sce <- scpModelWorkflow(
    sce,
    formula = ~ 1 + ## intercept
        ## normalisation
        MedianIntensity +
        ## batch effects
        PlateRow +
        Subject +
        ## biological variability
        HeartLocation,
    verbose = FALSE)

scpModelFilterThreshold(sce) <- 3

#########################################################
## Serialise object

save(ai2025a,
     file = file.path("../extdata/ai2025a.rda"),
     compress = "xz",
     compression_level = 9)