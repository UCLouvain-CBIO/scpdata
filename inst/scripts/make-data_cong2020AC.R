
####---- Cong et al. 2020, An. Chem. ----####


## Cong, Yongzheng, Yiran Liang, Khatereh Motamedchaboki, Romain 
## Huguet, Thy Truong, Rui Zhao, Yufeng Shen, Daniel Lopez-Ferrer, 
## Ying Zhu, and Ryan T. Kelly. 2020. “Improved Single-Cell Proteome 
## Coverage Using Narrow-Bore Packed NanoLC Columns and Ultrasensitive 
## Mass Spectrometry.” Analytical Chemistry, January. 
## https://doi.org/10.1021/acs.analchem.9b04631.

## The data files were donwloaded from: 
## ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2020/02/PXD016921

library(scp)
library(tidyverse)
setwd("inst/scripts/")
dataDir <- "../extdata/cong2020AC/"


####---- PSM data ----####


## Load the quantification data
list.files(path = dataDir,
           pattern = "evidence", 
           full.names = TRUE) %>%
  read.table(sep = "\t", header = TRUE) %>%
  ## Rename column
  dplyr::rename(Batch = Raw.file,  
                LFQ = Intensity) %>%
  ## Create 'Idenitfication.type' field as found in the peptide data 
  mutate(Identification.type = ifelse(Type == "MULTI-MATCH", 
                                      "By matching",
                                      "By MS/MS")) ->
  psms

## Create the sample metadata
meta <- data.frame(Batch = unique(psms$Batch), 
                   Channel = "LFQ",
                   CellNumber = c(20, 100, rep(1, 4), 0)) 

## Create the QFeatures object
cong2020AC <- readSCP(psms, 
                      meta, 
                      channelCol = "Channel", 
                      batchCol = "Batch")


####---- Peptide data ----####


## Load the quantification data
list.files(path = dataDir,
           pattern = "peptide", 
           full.names = TRUE) %>%
  read.table(sep = "\t", header = TRUE) ->
  pep

## Create the SingleCellExperiment object
pep <- readSingleCellExperiment(pep, 
                                ecol = grep("^Intensity[.]", 
                                            colnames(pep)), 
                                sep = "\t")

## Rename columns so they math with the PSM data
colnames(pep) %>%
  sub(pattern = "^Intensity[.]", replacement = "") %>%
  gsub(pattern = "[.]", replacement = " ") %>%
  paste0("_LFQ") ->
  colnames(pep)
## Name rows with peptide sequence
rownames(pep) <- rowData(pep)$Sequence

## Include the peptide data in the QFeatures object
cong2020AC <- addAssay(cong2020AC, pep, name = "peptides")

## Link the PSMs and the peptides
cong2020AC <- addAssayLink(cong2020AC, 
                           from = 1:7, 
                           to = "peptides", 
                           varFrom = rep("Sequence", 7),
                           varTo = "Sequence")


####---- Protein data ----####


## Load the quantification data
list.files(path = dataDir,
           pattern = "proteinGroups", 
           full.names = TRUE) %>%
  read.table(sep = "\t", header = TRUE) %>%
  ## Get a unique protein ID
  mutate(Protein = sub("^([^;]*);.*$", 
                       "\\1", 
                       Majority.protein.IDs)) ->
  prot

## Create the SingleCellExperiment object
prot <- readSingleCellExperiment(prot, 
                                 ecol = grep("^Intensity[.]", 
                                             colnames(prot)), 
                                 sep = "\t")
## Rename columns so they math with the PSM data
colnames(prot) %>%
  sub(pattern = "^Intensity[.]", replacement = "") %>%
  gsub(pattern = "[.]", replacement = " ") %>%
  paste0("_LFQ") ->
  colnames(prot)
## Name rows with peptide sequence
rownames(prot) <- rowData(prot)$Protein

## Include the protein data in the QFeatures object
cong2020AC <- addAssay(cong2020AC, prot, name = "proteins")

## Link the PSMs and the peptides
cong2020AC <- addAssayLink(cong2020AC, 
                           from = "peptides", 
                           to = "proteins", 
                           varFrom = "Leading.razor.protein",
                           varTo = "Protein")


# Save data as Rda file
# Note: saving is assumed to occur in "scpdata/inst/scripts"
save(cong2020AC, 
     file = file.path("../extdata/scpdata/cong2020AC.Rda"),
     compress = "xz", 
     compression_level = 9)
