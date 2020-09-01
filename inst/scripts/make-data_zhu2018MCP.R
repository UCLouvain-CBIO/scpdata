
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
     file = file.path("../EHdata/zhu2018MCP.rda"),
     compress = "xz", 
     compression_level = 9)



####---- Experiment metadata ----####

# Create the experiment data for the MSnSet object
expdat <- new("MIAPE",
              title = "Spatially Resolved Proteome Mapping of Laser Capture Microdissected Tissue with Automated Sample Transfer to Nanodroplets",
              abstract = "Current mass spectrometry (MS)-based proteomics approaches are ineffective for mapping protein expression in tissue sections with high spatial resolution because of the limited overall sensitivity of conventional workflows. Here we report an integrated and automated method to advance spatially resolved proteomics by seamlessly coupling laser capture microdissection (LCM) with a recently developed nanoliter-scale sample preparation system termed nanoPOTS (Nanodroplet Processing in One pot for Trace Samples). The workflow is enabled by prepopulating nanowells with DMSO, which serves as a sacrificial capture liquid for microdissected tissues. The DMSO droplets efficiently collect laser-pressure catapulted LCM tissues as small as 20 μm in diameter with success rates >87%. We also demonstrate that tissue treatment with DMSO can significantly improve proteome coverage, likely due to its ability to dissolve lipids from tissue and enhance protein extraction efficiency. The LCM-nanoPOTS platform was able to identify 180, 695, and 1827 protein groups on average from 12-μm-thick rat brain cortex tissue sections having diameters of 50, 100, and 200 μm, respectively. We also analyzed 100-μm-diameter sections corresponding to 10–18 cells from three different regions of rat brain and comparatively quantified ∼1000 proteins, demonstrating the potential utility for high-resolution spatially resolved mapping of protein expression in tissues.",
              url = "http://dx.doi.org/10.1074/mcp.TIR118.000686",
              dateStamp = "01/10/2018",
              name = "Zhu, Ying; Dou, Maowei; Piehowski, Paul D; Liang, Yiran; Wang, Fangjun; Chu, Rosalie K; Chrisler, William B; Smith, Jordan N; Schwarz, Kaitlynn C; Shen, Yufeng; Shukla, Anil K; Moore, Ronald J; Smith, Richard D; Qian, Wei-Jun; Kelly, Ryan T",
              lab = "Environmental Molecular Sciences Laboratory, Pacific Northwest National Laboratory",
              instrumentModel = "Orbitrap Fusion Lumos Tribrid mass spectrometer",
              instrumentManufacturer = "ThermoFisher",
              softwareName = "",
              softwareVersion = "",
              switchingCriteria = "",
              ionSource = "",
              ionSourceDetails = "",
              analyser = "",
              analyserDetails = "",
              collisionEnergy = "")


