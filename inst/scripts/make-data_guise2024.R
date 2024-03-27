####---- Guise et al. 2024 ---####


## Guise, Amanda J., Santosh A. Misal, Richard Carson, Jen-Hwa Chu,
## Hannah Boekweg, Daisha Van Der Watt, Nora C. Welsh, et al. 2024.
## “TDP-43-Stratified Single-Cell Proteomics of Postmortem Human
## Spinal Motor Neurons Reveals Protein Dynamics in Amyotrophic
## Lateral Sclerosis.” Cell Reports 43 (1): 113636.

## All files were downloaded from the MASSIVE repo MSV000092119

library(SingleCellExperiment)
library(scp)
library(tidyverse)
dataDir <-"~/Documents/.localData/SCP/guise2024/"

####---- Load sample annotations ----####

runInfo <- read.table(
    list.files(dataDir, "InputFiles", full.names = TRUE),
    header = TRUE, sep = "\t"
)
colAnnot <- read.table(
    list.files(dataDir, "^Groups.txt", full.names = TRUE),
    header = TRUE, sep = "\t"
)
colnames(colAnnot)[1] <- "File.ID"
colAnnot <- merge(
    colAnnot, runInfo[, c("File.ID", "File.Name", "Creation.Date")],
    by = "File.ID"
)
## Add the colum name that contains quantification columns in the PSM
## table
colAnnot$QuantCol <- "Precursor.Abundance"

####---- Load PSM data ----####

psm <- read.table(
    list.files(dataDir, "PSMs.txt", full.names = TRUE),
    header = TRUE, sep = "\t"
)
## Create QFeatures object
guise2024 <- readSCP(
    psm, colAnnot, batchCol = "File.ID", channelCol = "QuantCol",
    suffix = "", sep = ""
)
## Warning: Missing metadata. The features are removed for F61, F34,
## F42, F88, F77, F8, F21, F5

####---- Load peptide data ----####

pep <- read.table(
    list.files(dataDir, "PeptideG", full.names = TRUE),
    header = TRUE, sep = "\t"
)
pep <- readSingleCellExperiment(
    pep, ecol = grep("^Abundances..Grouped...F", colnames(pep)),
    fnames = "Peptide.Groups.Peptide.Group.ID"
)
## Adapt the column names to match those in the QFeatuers object
colnames(pep) <- sub("^Ab.*(F\\d+).*$", "\\1", colnames(pep))
## Add to QFeatures
guise2024 <- addAssay(guise2024, pep, "peptides")
## NOTE: could not map psms to peptides because I could not find a
## variable that uniquely maps the two levels.

####---- Load protein data ----####

prot <- read.table(
    list.files(dataDir, "Proteins.txt", full.names = TRUE),
    header = TRUE, sep = "\t"
)
prot <- readSingleCellExperiment(
    prot, ecol = grep("^Abundance..F", colnames(prot)),
    fnames = "Accession"
)
## Adapt the column names to match those in the QFeatuers object
colnames(prot) <- sub("^Ab.*(F\\d+).*$", "\\1", colnames(prot))
## Add to QFeatures
guise2024 <- addAssay(guise2024, prot, "proteins")
## NOTE: could not map psms to proteins because some psms have multiple
## associated proteins (`Master.Protein.Accessions`). This is not yet
## supported in `QFeatures`.

####---- Save dataset ----####


## Save the `QFeatures` object in an Rda file that will be uploaded to
## ExperimentHub
save(guise2024,
     compress = "xz",
     compression_level = 9,
     file = "~/Documents/.localData/scpdata/guise2024.rda")

