
####---- Zhu et al. 2020, eLife ----####


## Zhu, Ying, Mirko Scheibinger, Daniel Christian Ellwanger, Jocelyn F. 
## Krey, Dongseok Choi, Ryan T. Kelly, Stefan Heller, and Peter G. 
## Barr-Gillespie. 2019. “Single-Cell Proteomics Reveals Changes in 
## Expression during Hair-Cell Development.” eLife 8 (November). 
## https://doi.org/10.7554/eLife.50777.

## The data files were donwloaded from: 
## ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/11/PXD014256
## 
## PSMS, peptide and proteins data are found in
## SEARCH.zip -> SEARCH/Experiment 1 + 2/txt 2019-05-09a/
## Sample annotation is found in 
## OTHER.zip -> OTHER/Zhu_2019_chick_single_cell_samples.xlsx

library(openxlsx)
library(scp)
library(tidyverse)
setwd("inst/scripts/")
dataDir <- "../extdata/zhu2019EL/"


####---- PSM data ----####


list.files(path = dataDir,
           pattern = "samples.xlsx",
           full.names = TRUE) %>%
  ## The annotation is in the 3rd sheet
  read.xlsx(sheet = 3, colNames = TRUE, startRow = 7) %>%
  ## Rename column to match with PSM data
  rename(Raw.file = RAW.file.name) %>%
  ## Add channel that will point to the quantification column
  mutate(Channel = "LFQ") -> 
  meta
  

####---- PSM data ----####


## Load the quantification data
list.files(path = dataDir,
           pattern = "evidence", 
           full.names = TRUE) %>%
  read.table(sep = "\t", header = TRUE) %>%
  ## Rename column to match the annotation table
  rename(LFQ = Intensity) %>%
  ## Add file extension to file name to match metadata
  mutate(Raw.file = paste0(Raw.file, ".raw")) %>%
  ## Select only batches that are annotated
  filter(Raw.file %in% meta$Raw.file) %>%
  ## Add the sample name
  left_join(meta[, c("Raw.file", "Sample.name")], 
            by = "Raw.file") ->
  psms

## Create the QFeatures object
zhu2019EL <- readSCP(psms, 
                     meta, 
                     channelCol = "Channel", 
                     batchCol = "Sample.name")


####---- Peptide data ----####


## Load the quantification data
list.files(path = dataDir,
           pattern = "peptide", 
           full.names = TRUE) %>%
  read.table(sep = "\t", header = TRUE) ->
  pep

## Rename columns so they math with the PSM data
colnames(pep) %>%
  sub(pattern = "^Intensity.(.*)$", replacement = "\\1_LFQ") ->
  colnames(pep)

## Create the SingleCellExperiment object
pep <- readSingleCellExperiment(pep, 
                                ecol = grep("LFQ$", colnames(pep)))

## Name rows with peptide sequence
rownames(pep) <- rowData(pep)$Sequence

## Include the peptide data in the QFeatures object
zhu2019EL <- addAssay(zhu2019EL, pep, name = "peptides")
## Link the PSMs and the peptides
zhu2019EL <- addAssayLink(zhu2019EL, 
                           from = 1:60, 
                           to = "peptides", 
                           varFrom = rep("Sequence", 60),
                           varTo = "Sequence")


####---- Protein data ----####


## Load the quantification data
list.files(path = dataDir,
           pattern = "proteinGroups.txt",
           full.names = TRUE) %>%
  read.table(sep = "\t", header = TRUE) %>%
  ## Get a unique protein ID
  mutate(Protein = sub("^([^;]*);.*$", 
                       "\\1", 
                       Majority.protein.IDs)) ->
  prot

## Rename columns so they math with the PSM data
colnames(prot) %>%
  sub(pattern = "^Intensity.(.*)$", replacement = "\\1_LFQ") ->
  colnames(prot)


## Create the SingleCellExperiment object
prot <- readSingleCellExperiment(prot, 
                                 ecol = grep("LFQ$", colnames(prot)))
## Name rows with peptide sequence
rownames(prot) <- rowData(prot)$Protein

## Include the protein data in the QFeatures object
zhu2019EL <- addAssay(zhu2019EL, prot, name = "proteins")
## Link the PSMs and the peptides
zhu2019EL <- addAssayLink(zhu2019EL, 
                          from = "peptides", 
                          to = "proteins", 
                          varFrom = "Leading.razor.protein",
                          varTo = "Protein")


# Save data as Rda file
# Note: saving is assumed to occur in "scpdata/inst/scripts"
save(zhu2019EL, 
     file = file.path("../EHdata/zhu2019EL.Rda"),
     compress = "xz", 
     compression_level = 9)

