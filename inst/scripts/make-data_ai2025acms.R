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

save(ai2025a,
     file = file.path("../extdata/ai2025a.rda"),
     compress = "xz",
     compression_level = 9)