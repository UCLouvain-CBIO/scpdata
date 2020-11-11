
####---- Zhu et al. 2018, Nature Communications - HeLa dilutions ----####

## Zhu, Ying, Paul D. Piehowski, Rui Zhao, Jing Chen, Yufeng Shen, Ronald J. 
## Moore, Anil K. Shukla, et al. 2018. “Nanodroplet Processing Platform for Deep 
## and Quantitative Proteome Profiling of 10-100 Mammalian Cells.” Nature 
## Communications 9 (1): 882.

library(SingleCellExperiment)
library(scp)
library(tidyverse)
setwd("./inst/scripts")
dataDir <- "../extdata/zhu2018NC_hela/"

## The data was downloaded from ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2018/01/PXD006847
## to scpdata/inst/extdata/zhu2018NC_hela


####---- Peptide data ----####


## Load the quantification data
paste0(dataDir, "CulturedCells_peptides.txt") %>%
  read.table(sep = "\t", header = TRUE) %>%
  mutate(Batch = "peptides") ->
  dat

## Create the sample metadata
samples <- grep("Intensity[.]", colnames(dat), value = TRUE)
data.frame(Batch = "peptides",
           Sample = samples) %>%
  mutate(SampleType = sub("Intensity[.]", "", Sample),
         SampleType = ifelse(grepl("blank", SampleType), "Blank", SampleType),
         SampleType = ifelse(grepl("Lysate", SampleType), "Lysate", SampleType),
         SampleType = ifelse(grepl("cel", SampleType), "Hela", SampleType),
         SampleType = ifelse(grepl("MCF7", SampleType), "MCF7", SampleType),
         SampleType = ifelse(grepl("THP1", SampleType), "THP1", SampleType),
         CellNumber = ifelse(SampleType == "Hela", 
                             sub("^.*Hela(.*)cel.*$", "\\1", Sample), 
                             NA),
         CellNumber = ifelse(SampleType == "Lysate", 50, CellNumber),
         CellNumber = as.numeric(CellNumber)) ->
  meta

## Create the QFeatures object
zhu2018NC_hela <- readSCP(dat, 
                          meta, 
                          channelCol = "Sample", 
                          batchCol = "Batch")

# Save data as Rda file
# Note: saving is assumed to occur in "(...)/scpdata/inst/scripts"
save(zhu2018NC_hela, 
     file = file.path("../extdata/scpdata/zhu2018NC_hela.Rda"),
     compress = "xz", 
     compression_level = 9)




