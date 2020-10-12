
####---- Zhu et al. 2018, Molecular & Cellular Proteomics ----####

# Zhu, Ying, Maowei Dou, Paul D. Piehowski, Yiran Liang, Fangjun Wang, Rosalie 
# K. Chu, William B. Chrisler, et al. 2018. “Spatially Resolved Proteome Mapping 
# of Laser Capture Microdissected Tissue with Automated Sample Transfer to 
# Nanodroplets.” Molecular & Cellular Proteomics: MCP 17 (9): 1864–74.

library(scp)
library(mzR)
library(tidyverse)
setwd("inst/scripts/")
dataDir <- "../extdata/zhu2018MCP/"

# The data was downloaded from ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2018/07/PXD008844
# to scpdata/inst/extdata/zhu2018MCP


####---- Peptide data ----####


## Load the peptide data
list.files(path = dataDir,
           pattern = "Peptides",
           full.names = TRUE) %>%
  read.table(sep = "\t", header = TRUE) %>%
  mutate(Batch = "peptides") ->
  quant

## Create the metadata table 
data.frame(Batch = "peptides",
           Channel = grep("^Intensity.", colnames(quant), value = TRUE)) %>%
  mutate(SectionWidth = sub("^.*[.]([0-9]*um)_.*$", "\\1", Channel),
         SampleType = str_extract(Channel, "CC|CP|Mix"),
         SampleType = ifelse(is.na(SampleType), "CTX", SampleType),
         Prep = str_match(Channel, "08.*Prep$"),
         Replicate = str_match(Channel, "(_|Mix)(\\d)(_|$)")[, 3]) ->
  meta 

## Create the QFeatures object
zhu2018MCP <- readSCP(quant,
                      meta, 
                      batchCol = "Batch", 
                      channelCol = "Channel")

### Save data as Rda file
## Note: saving is assumed to occur in "scpdata/inst/scripts"
save(zhu2018MCP, 
     file = file.path("../EHdata/scpdata/zhu2018MCP.Rda"),
     compress = "xz", 
     compression_level = 9)
