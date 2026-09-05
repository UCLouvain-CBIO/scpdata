library("SCP.replication")
library("tidyverse")
library("data.table")

## report <- fread("~/PhD/.localdata/SCP/brunner2022/20210919_DIANN_SingleCellOutput.tsv")
report <- fread("~/Downloads/DIANN1.8_SingleCells_CellCycle/20210919_DIANN_SingleCellOutput.tsv")
report$File.Name <- sub("^.+202010", "202010", report$File.Name)

## Build the annotation table from the run names
annot <- DataFrame(File.Name = unique(report$File.Name))
annot$Run <- gsub("^D.*[\\]|[.]d$", "", annot$File.Name)
annot$runCol <- annot$File.Name
otherVars <- strsplit(annot$Run, "_")
otherVars <- lapply(otherVars, function(x) {
    if (length(x) == 12) x <- x[-4]
    x
})
otherVars <- do.call(rbind, otherVars)
colnames(otherVars) <- c("Date", "MsInstrument", "Purification", "User",
                         "SampleAnnotation", "SampleType", "CellCycleStage",
                         "..undetermined..", "PlatePosition", "CellNumber",
                         "RunID")
annot <- cbind(annot, otherVars)

## Format to a QFeatures object
## brunner2022 <- readSCPfromDIANN(annot, report)
brunner2022 <- readSCPfromDIANN(assayData = report, colData = annot)
names(brunner2022) <- sub("^.+2020", "2020", names(brunner2022))

brunner2022 <- aggregateFeatures(
    brunner2022,
    i = seq_along(brunner2022),
    fcol = "Modified.Sequence",
    name = paste0("modpep_",
                  sub("\\.d", "", sub("^.+_", "", names(brunner2022)))),
    fun = colMedians,
    na.rm = TRUE)

x <- joinAssays(brunner2022,
                grep("modpep", names(brunner2022)),
                name = "modpep")

## The protein data
## pgTable <- read.delim("~/PhD/.localdata/SCP/brunner2022/20210919_DIANN_SingleCellOutput.pg_matrix.tsv",
##                       check.names = FALSE)

pgTable <- read.delim("~/Downloads/DIANN1.8_SingleCells_CellCycle/20210919_DIANN_SingleCellOutput.pg_matrix.tsv",
                      check.names = FALSE)
prots <- readSingleCellExperiment(pgTable, quantCols = unique(report$File.Name),
                                  fnames = "Protein.Names")
prots <- prots[rownames(prots) != ""]
rownames(prots) <- make.unique(rownames(prots))
from <- names(brunner2022)
brunner2022 <- addAssay(brunner2022, prots, "proteins")
brunner2022 <- addAssayLink(brunner2022,
                            from = from,
                            to = "proteins",
                            varFrom = rep("Protein.Group", length(from)),
                            varTo = "Protein.Group")

# Save data as Rda file
save(brunner2022,
     file = "~/PhD/.localdata/scpdata/brunner2022.Rda",
     compress = "xz",
     compression_level = 9)
