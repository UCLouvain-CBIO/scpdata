
####---- Zhu et al. 2018, Nature Communications - HeLa lysates ----####

## Zhu, Ying, Paul D. Piehowski, Rui Zhao, Jing Chen, Yufeng Shen, Ronald J. 
## Moore, Anil K. Shukla, et al. 2018. “Nanodroplet Processing Platform for Deep 
## and Quantitative Proteome Profiling of 10-100 Mammalian Cells.” Nature 
## Communications 9 (1): 882.

library(SingleCellExperiment)
library(scp)
library(tidyverse)
setwd("./inst/scripts")
dataDir <- "../extdata/zhu2018NC_lysates/"

## The data was downloaded from ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2018/01/PXD006847
## to scpdata/inst/extdata/zhu2018NC_lysates


####---- Peptide data ----####


## Load the quantification data
paste0(dataDir, "Vail_Prep_Vail_peptides.txt") %>%
  read.table(sep = "\t", header = TRUE) %>%
  mutate(Batch = "peptides") ->
  dat

## Create the sample metadata
samples <- grep("Intensity[.]", colnames(dat), value = TRUE)
data.frame(Batch = "peptides",
           Column = samples) %>%
  mutate(Sample = sub("Intensity[.]", "", Column),
         CellEquivalent = sub("^([0-9]*)HeLa.*$", "\\1", Sample),
         Replicate = sub("^.*([0-9])$", "\\1", Sample),
         Digestion = sub("^.*HeLa_(.*)_08.*$", "\\1", Sample),
         Date = sub("^.*p_(.*)_.$", "\\1", Sample),
         Date = as.Date(Date, "%d%m%y")) ->
  meta

## Create the QFeatures object
zhu2018NC_lysates <- readSCP(dat, 
                             meta, 
                             channelCol = "Column", 
                             batchCol = "Batch")

# Save data as Rda file
# Note: saving is assumed to occur in "(...)/scpdata/inst/scripts"
save(zhu2018NC_lysates, 
     file = file.path("../EHdata/scpdata/zhu2018NC_lysates.Rda"),
     compress = "xz", 
     compression_level = 9)




