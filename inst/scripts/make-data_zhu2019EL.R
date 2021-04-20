
####---- Zhu et al. 2020, eLife ----####


## Zhu, Ying, Mirko Scheibinger, Daniel Christian Ellwanger, Jocelyn F. 
## Krey, Dongseok Choi, Ryan T. Kelly, Stefan Heller, and Peter G. 
## Barr-Gillespie. 2019. “Single-Cell Proteomics Reveals Changes in 
## Expression during Hair-Cell Development.” eLife 8 (November). 
## https://doi.org/10.7554/eLife.50777.

## The data files were downloaded from: 
## ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/11/PXD014256
## 
## PSMS, peptide and proteins data are found in
## SEARCH.zip -> SEARCH/Experiment 1 + 2/txt 2019-05-09a/
## Sample annotation is found in 
## OTHER.zip -> OTHER/Zhu_2019_chick_single_cell_samples.xlsx

library(openxlsx)
library(scp)
library(tidyverse)
dataDir <- "../.localdata/SCP/zhu2019EL/"


####---- Sample annotation ----####


list.files(path = dataDir,
           pattern = "samples_CORRECTED.xlsx",
           full.names = TRUE) %>%
    ## The annotation is in the 3rd sheet
    read.xlsx(sheet = 3, colNames = TRUE, startRow = 7) %>%
    ## Rename column to match with PSM data
    rename(Raw.file = RAW.file.name,
           ## Avoid issue with special characters
           FM1.43.signal = `FM1-43.signal`) %>%
    ## Add the column name containing the quantitative data
    mutate(QuantCol = "Intensity",
           ## Add the experimental replicate
           Experiment = sub("^(.).*$", "\\1", Sample.name)) -> 
    meta


####---- PSM data ----####


## Load the quantification data
list.files(path = dataDir,
           pattern = "evidence", 
           full.names = TRUE) %>%
    read.table(sep = "\t", header = TRUE) %>%
    ## Add file extension to file name to match metadata
    mutate(Raw.file = paste0(Raw.file, ".raw")) %>% 
    ## Select only batches that are annotated
    filter(Raw.file %in% meta$Raw.file) %>%
    ## Add the sample name
    left_join(meta[, c("Raw.file", "Sample.name")], 
              by = "Raw.file") ->
    psms
## Note that we miss annotations for the following runs: 
## - Single_Hair_Cell_OHSU_1cell_Low_030819_R10_YF30um_350bar.raw
## - Single_Hair_Cell_OHSU_1cell_High_030819_R10_YF30um_350bar.raw

## Create the QFeatures object
zhu2019EL <- readSCP(featureData = psms, 
                     colData = meta, 
                     channelCol = "QuantCol", 
                     batchCol = "Sample.name",
                     suffix = "")


####---- Peptide data ----####


## Load the quantification data
list.files(path = dataDir,
           pattern = "peptide", 
           full.names = TRUE) %>%
    read.table(sep = "\t", header = TRUE) ->
    pep

## Rename columns so they match with the PSM data
colnames(pep) <- sub(pattern = "^Intensity.", 
                     replacement = "",
                     colnames(pep))

## Create the SingleCellExperiment object
pep <- readSingleCellExperiment(pep, 
                                ecol = meta$Sample.name)

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

## Rename columns so they match with the PSM data
colnames(prot) <- sub(pattern = "^Intensity.", 
                      replacement = "",
                      colnames(prot))

## Create the SingleCellExperiment object
prot <- readSingleCellExperiment(prot, 
                                 ecol = meta$Sample.name)
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
     file = file.path("../.localdata/scpdata/zhu2019EL.Rda"),
     compress = "xz", 
     compression_level = 9)

