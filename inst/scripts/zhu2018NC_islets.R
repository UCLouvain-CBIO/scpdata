
####---- Zhu et al. 2018, Nature Communications - T1D islets ----####

## Zhu, Ying, Paul D. Piehowski, Rui Zhao, Jing Chen, Yufeng Shen, Ronald J. 
## Moore, Anil K. Shukla, et al. 2018. “Nanodroplet Processing Platform for Deep 
## and Quantitative Proteome Profiling of 10-100 Mammalian Cells.” Nature 
## Communications 9 (1): 882.

library(MSnbase)
library(mzR)
library(SingleCellExperiment)
setwd("inst/scripts/")

# The data was downloaded from ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2018/01/PXD006847
# to scpdata/extdata/zhu2018NC


####---- Experiment metadata ----####

# Create the experiment data for the MSnSet object
expdat <- new("MIAPE",
              title = "Nanodroplet Processing Platform for Deep and Quantitative Proteome Profiling of 10-100 Mammalian Cells.",
              abstract = "Nanoscale or single-cell technologies are critical for biomedical applications. However, current mass spectrometry (MS)-based proteomic approaches require samples comprising a minimum of thousands of cells to provide in-depth profiling. Here, we report the development of a nanoPOTS (nanodroplet processing in one pot for trace samples) platform for small cell population proteomics analysis. NanoPOTS enhances the efficiency and recovery of sample processing by downscaling processing volumes to <200 nL to minimize surface losses. When combined with ultrasensitive liquid chromatography-MS, nanoPOTS allows identification of ~1500 to ~3000 proteins from ~10 to ~140 cells, respectively. By incorporating the Match Between Runs algorithm of MaxQuant, >3000 proteins are consistently identified from as few as 10 cells. Furthermore, we demonstrate quantification of ~2400 proteins from single human pancreatic islet thin sections from type 1 diabetic and control donors, illustrating the application of nanoPOTS for spatially resolved proteome measurements from clinical tissues.",
              url = "http://dx.doi.org/10.1038/s41467-018-03367-w",
              dateStamp = "05/07/2017",
              name = "Zhu, Ying, Paul D. Piehowski, Rui Zhao, Jing Chen, Yufeng Shen, Ronald J. Moore, Anil K. Shukla, et al",
              lab = "PNNL",
              instrumentModel = "Orbitrap Fusion Lumos Tribrid MS ",
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
f <- "../../extdata/zhu2018NC/Islet_t1d_ct_peptides.txt"
colSel <- grepEcols(f, pattern = "Intensity.", split = "\t")
.fnames <- c("Sequence", "Mass", "Proteins", "Leading.razor.protein", 
             "Gene.names", "Protein.names", "Charges", "PEP", "Score", 
             "Reverse", "Potential.contaminant", "Missed cleavages")
dat <- readMSnSet2(f, sep = "\t", ecol = colSel, fnames = .fnames)

## Use columns for generating phenotype data
pData(dat) <- data.frame(row.names = sampleNames(dat))
pData(dat)$sampleType <- "control"
pData(dat)$sampleType[grepl("[.]T", sampleNames(dat))] <- "T1D"

# Remove contaminants and reverse hits
dat <- dat[fData(dat)$Reverse != "+"]
dat <- dat[fData(dat)$Potential.contaminant != "+"]
# Remove bad identification
dat <- dat[fData(dat)$PEP < 0.01]

# Replace 0 by NAs
ed <- exprs(dat)
ed[ed == 0] <- NA
exprs(dat) <- ed

# Create the SingleCellExperiment object
zhu2018NC_islets_peptide <-
  SingleCellExperiment(assay = list(peptide = exprs(dat)), 
                       colData = pData(dat),
                       rowData = fData(dat),
                       metadata = list(experimentData = expdat))
# Save data as Rda file
# Note: saving is assumed to occur in "(...)/scpdata/inst/scripts"
save(zhu2018NC_islets_peptide, compress = "xz", compression_level = 9,
     file = file.path("../../data/zhu2018NC_islets_peptide.rda"))


