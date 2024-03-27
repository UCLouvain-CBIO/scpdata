## Core packages of this workflow
library("scp")
library("SCP.replication")
library("tidyverse")

## All data files were retrieved from this Google Drive repository:
## https://drive.google.com/drive/folders/1pUC2zgXKtKYn22mlor0lmUDK0frgwL_-
datadir <- "~/PhD/.localdata/SCP/derks2022/"

## Sample annotations


sampleAnnot <- read.delim(paste0(datadir, "Meta_SingleCell_updated_1.tsv"))
## Add which dataset each sample is part of
sampleAnnot$dataset <- sampleAnnot$Instrument
sampleAnnot$dataset[!sampleAnnot$Real_single_cell & 
                        sampleAnnot$Instrument == "Q-Exactive"] <- "bulk"
## Adapt variable to better match the DIANN output data
sampleAnnot$Label <- as.character(sampleAnnot$Label)
timsPath <- make.names("F:\\JD\\plexDIA\\Bruker\\OneDrive_1_3-9-2022\\")
qePath <- make.names("F:\\JD\\plexDIA\\nPOP\\")
sampleAnnot$File.Name <- ifelse(sampleAnnot$Instrument == "timsTOFSCP",
                                paste0(timsPath, sampleAnnot$Raw, ".d"),
                                paste0(qePath, sampleAnnot$Raw, ".raw"))

## Bulk data

## cf https://drive.google.com/drive/folders/1yRzuIXnMbt-_8_skOgiVkJksqc20RoLK

# We load the DIA-NN main output table and the MS1 extracted report table.
# These are read and combined in a `QFeatures`object. 
extractedDataBulk <- read.delim(paste0(datadir, "qe_bulk/Report.pr_matrix_channels_ms1_extracted.tsv"))
reportDataBulk <- read.delim(paste0(datadir, "qe_bulk/Report.tsv"))
reportDataBulk$File.Name <- make.names(reportDataBulk$File.Name)
bulk <- readSCPfromDIANN(colData = sampleAnnot, 
                         reportData = reportDataBulk,
                         extractedData = extractedDataBulk,
                         multiplexing = "mTRAQ")
## Rename the MS1Extracted assay
names(bulk)[length(bulk)] <- "bulk_prec_extracted"

## Load timsTOF-SCP data

## cf https://drive.google.com/drive/folders/1RosRkMdYfhbQ-XtUNKly0TdCrZrUQVmO

# We load the DIA-NN main output table and the MS1 extracted report table.
# These are read and combined in a `QFeatures`object. 
extractedDataTims <- read.delim(paste0(datadir, "tims_sc/Report.pr_matrix_channels_ms1_extracted.tsv"))
reportDataTims <- read.delim(paste0(datadir, "tims_sc/Report.tsv"))
reportDataTims$File.Name <- make.names(reportDataTims$File.Name)
# We modify the `Run` variable to match the `Run` variables in the other tables
reportDataTims$Run <- make.names(reportDataTims$Run)
tims <- readSCPfromDIANN(colData = sampleAnnot, 
                         reportData = reportDataTims,
                         extractedData = extractedDataTims,
                         multiplexing = "mTRAQ")
names(tims)[length(tims)] <- "tims_prec_extracted"

## Load Q-Exactive data

## cf https://drive.google.com/drive/folders/1b_pavZ2sufR3oYBy2z9fVepCRDKXGvyF

# We load the DIA-NN main output table and the MS1 extracted report table.
# These are read and combined in a `QFeatures`object. 
extractedDataQE <- read.delim(paste0(datadir, "qe_sc/Report.pr_matrix_channels_ms1_extracted.tsv"))
reportDataQE <- read.delim(paste0(datadir, "qe_sc/Report.tsv"))
reportDataQE$File.Name <- make.names(reportDataQE$File.Name)
qe <- readSCPfromDIANN(colData = sampleAnnot, 
                       reportData = reportDataQE,
                       extractedData = extractedDataQE,
                       multiplexing = "mTRAQ")
names(qe)[length(qe)] <- "qe_prec_extracted"

## Load protein data

prots <- read.delim(paste0(datadir, "Proteins_SC_IDs.txt"))
prots <- readSingleCellExperiment(prots, fname = "prot",
                                  ecol = grep("id", colnames(prots)))
colData(prots) <- DataFrame(sampleAnnot[sampleAnnot$id %in% colnames(prots), ])
colnames(prots) <- paste0(prots$File.Name, ".", prots$Label)

## Combine all datasets
derks2022 <- c(bulk, tims, qe)
derks2022 <- addAssay(derks2022, prots, name = "proteins")
derks2022 <- addAssayLink(derks2022, 
                          from = grep("extracted$", names(derks2022)),
                          to = "proteins",
                          varFrom = rep("Protein.Group", 3), 
                          varTo = "prot")

## Save data
save(derks2022, 
     file = "~/PhD/.localdata/scpdata/derks2022.Rda",
     compress = "xz", 
     compression_level = 9)
