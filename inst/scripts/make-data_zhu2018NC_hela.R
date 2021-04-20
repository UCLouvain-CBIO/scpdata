
####---- Zhu et al. 2018, Nature Communications - HeLa dilutions ----####

## Zhu, Ying, Paul D. Piehowski, Rui Zhao, Jing Chen, Yufeng Shen, Ronald J. 
## Moore, Anil K. Shukla, et al. 2018. “Nanodroplet Processing Platform for Deep 
## and Quantitative Proteome Profiling of 10-100 Mammalian Cells.” Nature 
## Communications 9 (1): 882.

library(SingleCellExperiment)
library(scp)
library(tidyverse)
dataDir <- "~/PhD/.localdata/SCP/zhu2018NC_hela/"

## The data was downloaded from ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2018/01/PXD006847
## to scpdata/inst/extdata/zhu2018NC_hela

####---- Get the annotation data ----####

## The annotation file was manually created from the samples names and
## the available information from the methods section
meta <- read.csv(paste0(dataDir, "sample_annotation.csv"))
rownames(meta) <- paste0(meta$SampleName, "_sample")

####---- Peptide data ----####

## Load the quantification data
pep <-  read.table(paste0(dataDir, "CulturedCells_peptides.txt"),
                   sep = "\t", header = TRUE) 
colnames(pep) <- sub("Intensity.(.*)$", "\\1_sample", colnames(pep))
rownames(pep) <- pep$Sequence
pep <- readSingleCellExperiment(pep, ecol = grepl("sample$", colnames(pep)))

####---- Protein data ----####

## Load the quantification data
prot <-  read.table(paste0(dataDir, "CulturedCells_proteinGroups.txt"),
                    sep = "\t", header = TRUE)
colnames(prot) <- sub("Intensity.(.*)$", "\\1_sample", colnames(prot))
rownames(prot) <- prot$Protein.IDs
prot <- readSingleCellExperiment(prot, ecol = grepl("sample$", colnames(prot)))

## Create the QFeatures object
zhu2018NC_hela <- QFeatures(list(peptides = pep, proteins = prot), 
                            colData = DataFrame(meta))

# Save data as Rda file
# Note: saving is assumed to occur in "(...)/scpdata/inst/scripts"
save(zhu2018NC_hela, 
     file = "~/PhD/.localdata/scpdata/zhu2018NC_hela.Rda",
     compress = "xz", 
     compression_level = 9)
