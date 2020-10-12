
####---- Zhu et al. 2018, Nature Communications - T1D islets ----####

## Zhu, Ying, Paul D. Piehowski, Rui Zhao, Jing Chen, Yufeng Shen, Ronald J. 
## Moore, Anil K. Shukla, et al. 2018. “Nanodroplet Processing Platform for Deep 
## and Quantitative Proteome Profiling of 10-100 Mammalian Cells.” Nature 
## Communications 9 (1): 882.

library(SingleCellExperiment)
library(scp)
library(tidyverse)
setwd("./inst/scripts")
dataDir <- "../extdata/zhu2018NC_islets/"

# The data was downloaded from ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2018/01/PXD006847
# to scpdata/inst/extdata/zhu2018NC_islets


####---- Peptide data ----####


## Load the quantification data
paste0(dataDir, "Islet_t1d_ct_peptides.txt") %>%
  read.table(sep = "\t", header = TRUE) %>%
  mutate(Batch = "peptides") ->
  dat

## Create the sample metadata
samples <- grep("Intensity[.]", colnames(dat), value = TRUE)
data.frame(Batch = "peptides",
           Sample = samples) %>%
  mutate(SampleType = sub("^.*([CT]).*$", "\\1", Sample),
         SampleType = dplyr::recode(SampleType, 
                                    `T` = "T1D",
                                    `C` = "Control")) ->
  meta

## Create the QFeatures object
zhu2018NC_islets <- readSCP(dat, 
                            meta, 
                            channelCol = "Sample", 
                            batchCol = "Batch")

# Save data as Rda file
# Note: saving is assumed to occur in "(...)/scpdata/inst/scripts"
save(zhu2018NC_islets, 
     file = file.path("../EHdata/scpdata/zhu2018NC_islets.Rda"),
     compress = "xz", 
     compression_level = 9)
