
####---- Liang et al. 2020 ---####


## Liang, Yiran, Hayden Acor, Michaela A. McCown, Andikan J. Nwosu, 
## Hannah Boekweg, Nathaniel B. Axtell, Thy Truong, Yongzheng Cong, 
## Samuel H. Payne, and Ryan T. Kelly. 2020. “Fully Automated Sample 
## Processing and Analysis Workflow for Low-Input Proteome Profiling.”
##  Analytical Chemistry, December. 
##  https://doi.org/10.1021/acs.analchem.0c04240.

library(SingleCellExperiment)
library(scp)
library(tidyverse)
dataDir <- "~/PhD/.localdata/SCP/liang2020/"

## All data files were retrieved from the FTP links provided in the PRIDE 
## repository PXD021882
## https://www.ebi.ac.uk/pride/archive/projects/PXD021882

####---- Load sample annotation ----####

## Sample annotations were collected in a table from the methods 
## section and from table S3 in the paper
meta <- read.csv(paste0(dataDir, "HeLa_annotations.csv"))

####---- Load PSM data ----####

psms <- read.table(paste0(dataDir, "evidence.txt"), 
                   sep = "\t", 
                   header = TRUE)
    

####---- Create the QFeatures object ----####

liang2020_hela <- readSCP(featureData = psms, 
                             colData = meta, 
                             channelCol = "QuantificationColumn", 
                             batchCol = "Raw.file", 
                             suffix = "_sample")

####---- Load peptide data ----####

## Load the quantification data
pep <- read.table(paste0(dataDir, "peptides.txt"),
                  sep = "\t", 
                  header = TRUE)

## Rename columns so they match with the PSM data
colnames(pep) <- sub(pattern = "^Intensity.(.*)$", 
                     replacement = "HeLa_\\1_sample",
                     colnames(pep))
## Create the SingleCellExperiment object
pep <- readSingleCellExperiment(pep, 
                                ecol = grep("sample$", colnames(pep)))
## Name rows with peptide sequence
rownames(pep) <- rowData(pep)$Sequence
## Include the peptide data in the QFeatures object
liang2020_hela <- addAssay(liang2020_hela, 
                              pep, 
                              name = "peptides")
## Link the PSMs and the peptides
liang2020_hela <- addAssayLink(liang2020_hela, 
                                  from = 1:15, 
                                  to = "peptides", 
                                  varFrom = rep("Sequence", 15),
                                  varTo = "Sequence")

####---- Protein data ----####

## Load the quantification data
read.table(paste0(dataDir, "proteinGroups.txt"),
           sep = "\t",
           header = TRUE) %>% 
    ## Get a unique protein ID
    mutate(Protein = sub("^([^;]*);.*$", "\\1", 
                         Majority.protein.IDs)) ->
    prot

## Rename columns so they math with the PSM data
colnames(prot) <- sub(pattern = "^Intensity.(.*)$", 
                      replacement = "HeLa_\\1_sample", 
                      colnames(prot)) 

## Create the SingleCellExperiment object
prot <- readSingleCellExperiment(prot, 
                                 ecol = grep("sample$", colnames(prot)))
## Name rows with peptide sequence
rownames(prot) <- rowData(prot)$Protein

## Include the protein data in the QFeatures object
liang2020_hela <- addAssay(liang2020_hela, 
                              prot, 
                              name = "proteins")
## Link the PSMs and the peptides
liang2020_hela <- addAssayLink(liang2020_hela, 
                                  from = "peptides", 
                                  to = "proteins", 
                                  varFrom = "Leading.razor.protein",
                                  varTo = "Protein")

# Save data as Rda file
save(liang2020_hela, 
     file = "~/PhD/.localdata/scpdata/liang2020_hela.Rda",
     compress = "xz", 
     compression_level = 9)


