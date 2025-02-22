stopifnot(packageVersion("MsDataHub") >= "1.7.1")
stopifnot(packageVersion("QFeatures") >= "1.17.2")

library(tidyverse)
library(QFeatures)

ai2025iTab <- read_tsv(MsDataHub::Ai2025_iCMs_report.tsv())


tab <- tibble(File.Name = unique(ai2025iTab[[1]]))|>
    mutate(Sample = sub(".*\\\\", "",File.Name))|>
    mutate(Date = sub("^.*?(\\d{6}).*", "\\1", File.Name))|>
    mutate(Date = ymd(as.integer(Date)))|>
    mutate(Batch = sub("^.*?(Batch[0-9])_.+$", "\\1", File.Name))|>
    mutate(PlateWell = sub("^.*?([A-Z][0-9]+)_.+$", "\\1", File.Name))|>
    mutate(Position = sub("^.*?([A-Z]+[0-9]+)_1_[0-9]+\\.d$", "\\1", File.Name))|>
    mutate(Day = sub("^.*?(day[0-9]+)_.+$", "\\1", File.Name, ignore.case = T))|>
    mutate(Day = sub("\\D+", "\\1", Day))|>
    mutate(Day = as.numeric(Day))|>
    mutate(Cell_count = sub("^.*?([0-9]+cell).+$", "\\1", File.Name))

tab$runCol <- tab$File.Name

ai2025i <- readSCPfromDIANN(ai2025iTab
                            colData = DataFrame(tab),
                            fnames = "Precursor.Id")

## Setting the names of the QFeatures objects
ai2025i <- MultiAssayExperiment::renamePrimary(ai2025i, ai2025i$Sample)
names(ai2025i) <- ai2025i$Sample