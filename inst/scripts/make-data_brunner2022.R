library("tidyverse")
library("data.table")
library("QFeatures")
library("scp")

## The data files used below are available in
## DIANN1.8_SingleCells_CellCycle[.7z|zip] from PRIDE at
## https://www.ebi.ac.uk/pride/archive/projects/PXD024043?

report <- fread("~/tmp/DIANN1.8_SingleCells_CellCycle/20210919_DIANN_SingleCellOutput.tsv")
report$File.Name <- sub("^.+2020", "2020", report$File.Name) ## clean up filenames

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
## change set names
names(brunner2022) <- sub("\\.d", "", sub("2020.+_", "prec_", names(brunner2022)))

## Aggregate precursor into modified peptides
brunner2022 <-
    aggregateFeatures(brunner2022,
                      i = seq_along(brunner2022),
                      fcol = "Modified.Sequence",
                      name = sub("prec", "modpep", names(brunner2022)),
                      fun = colMedians,
                      na.rm = TRUE)

## Join modififed peptides into a wide table
brunner2022 <-
    joinAssays(brunner2022,
               grep("modpep", names(brunner2022)),
               name = "modpep")

## Log-transforme the modified peptide table
brunner2022 <-
    logTransform(brunner2022,
                 i = "modpep",
                 name = "logmodpep")

## Aggregate (log-transformed) modified peptides into stripped peptide sequences
brunner2022 <-
    aggregateFeatures(brunner2022,
                      i = "logmodpep",
                      name = "logpep",
                      fcol = "Stripped.Sequence",
                      fun = colMedians)
## Aggregate stripped peptides into (log-transformed) proteins
brunner2022 <-
    aggregateFeatures(brunner2022,
                      i = "logpep",
                      name = "logproteins",
                      fcol = "Protein.Group",
                      fun = colMedians)

## Add the protein data provided by the authors, and link it to the precursors
pgTable <- read.delim("~/tmp/DIANN1.8_SingleCells_CellCycle/20210919_DIANN_SingleCellOutput.pg_matrix.tsv",
                      check.names = FALSE)
names(pgTable) <- sub("^.+2020", "2020", names(pgTable))
prots <-
    readSingleCellExperiment(pgTable,
                             quantCols = unique(report$File.Name),
                             fnames = "Protein.Names")

prots <- prots[rownames(prots) != ""]
rownames(prots) <- make.unique(rownames(prots))
from <- grep("prec", names(brunner2022), value = TRUE)

brunner2022 <- addAssay(brunner2022, prots, "proteinsDiann")

brunner2022 <-
    addAssayLink(brunner2022,
                 from = from,
                 to = "proteinsDiann",
                 varFrom = rep("Protein.Group", length(from)),
                 varTo = "Protein.Group")

# Save data as Rda file
save(brunner2022,
     file = "~/tmp/DIANN1.8_SingleCells_CellCycle/brunner2022.Rda",
     compress = "xz",
     compression_level = 9)
