
####---- Zhu et al. 2018, Molecular & Cellular Proteomics ----####

# Zhu, Ying, Maowei Dou, Paul D. Piehowski, Yiran Liang, Fangjun Wang, Rosalie 
# K. Chu, William B. Chrisler, et al. 2018. “Spatially Resolved Proteome Mapping 
# of Laser Capture Microdissected Tissue with Automated Sample Transfer to 
# Nanodroplets.” Molecular & Cellular Proteomics: MCP 17 (9): 1864–74.

library(MSnbase)
library(SingleCellExperiment)
library(stringr)
setwd("inst/scripts/")

# The data was downloaded from ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2018/07/PXD008844
# to scpdata/extdata/zhu2018MCP


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


####---- Peptide data ----####

## Load the quantification data
f <- "../../extdata/zhu2018MCP/MaxQuant_Peptides.txt"
ecol <- grepEcols(f, pattern = "Intensity.", split = "\t")
fnames <- c("Sequence", "Length", "Mass", "Proteins", "Leading.razor.protein", 
             "Gene.names", "Protein.names", "Charges", "PEP", "Score", 
             "Reverse", "Potential.contaminant", "Missed.cleavages")
dat <- readMSnSet2(f, sep = "\t", ecol = ecol, fnames = fnames)

## Use columns for generating phenotype data
pData(dat) <- data.frame(row.names = sampleNames(dat))
pData(dat)$sectionWidth <- str_match(sampleNames(dat), "[.](.*)um")[, 2]
pData(dat)$sectionWidth <- as.numeric(pData(dat)$sectionWidth)
pData(dat)$sampleType <- str_extract(sampleNames(dat), "CC|CP|Mix")
pData(dat)$sampleType[is.na(pData(dat)$sampleType)] <- "CTX"
pData(dat)$Prep <- str_match(sampleNames(dat), "08.*Prep$")
pData(dat)$replicate <- str_match(sampleNames(dat), "(_|Mix)(\\d)(_|$)")[, 3]

# Remove contaminants and reverse hits
dat <- dat[fData(dat)$Reverse != "+"]
dat <- dat[fData(dat)$Potential.contaminant != "+"]
# Remove bad identification
dat <- dat[fData(dat)$PEP < 0.01]

# Replace 0 by NAs
ed <- exprs(dat)
ed[ed == 0] <- NA
exprs(dat) <- ed

# Remove rows with all NA's
dat <- filterNA(dat, pNA = 1-1E-6)

# Create the SingleCellExperiment object
zhu2018MCP_peptide <-
  SingleCellExperiment(assay = list(peptide = exprs(dat)), 
                       colData = pData(dat),
                       rowData = fData(dat),
                       metadata = list(experimentData = expdat))
# Save data as Rda file
# Note: saving is assumed to occur in "(...)/scpdata/inst/scripts"
save(zhu2018MCP_peptide, compress = "xz", compression_level = 9,
     file = file.path("../../data/zhu2018MCP_peptide.rda"))
