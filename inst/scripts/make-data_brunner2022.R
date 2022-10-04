library("SCP.replication")
library("tidyverse")
library("data.table")

report <- fread("~/PhD/.localdata/SCP/brunner2022/20210919_DIANN_SingleCellOutput.tsv")

## Build the annotation table from the run names
annot <- DataFrame(File.Name = unique(report$File.Name))
annot$Run <- gsub("^D.*[\\]|[.]d$", "", annot$File.Name)
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
brunner2022 <- readSCPfromDIANN(annot, report)

## The protein data 
pgTable <- read.delim("~/PhD/.localdata/SCP/brunner2022/20210919_DIANN_SingleCellOutput.pg_matrix.tsv", 
                      check.names = FALSE)
prots <- readSingleCellExperiment(pgTable, ecol = unique(report$File.Name), 
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
