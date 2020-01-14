
####---- Specht et al. 2019 ---####


# Specht, Harrison, Edward Emmott, Toni Koller, and Nikolai Slavov. 2019. 
# “High-Throughput Single-Cell Proteomics Quantifies the Emergence of Macrophage 
# Heterogeneity.” bioRxiv. https://doi.org/10.1101/665307.

library(MSnbase)
library(tidyr)
library(SingleCellExperiment)
setwd("inst/scripts")

# The data was downloaded from https://scope2.slavovlab.net/docs/data to 
# scpdata/inst/extdata/specht2019


####---- Experiment metadata ----####

# Create the experiment data for the MSnSet object
expdat <- new("MIAPE",
              title = "High-throughput single-cell proteomics quantifies the emergence of macrophage heterogeneity",
              abstract = "The fate and physiology of individual cells are controlled by networks of proteins. Yet, our ability to quantitatively analyze protein networks in single cells has remained limited. To overcome this barrier, we developed SCoPE2. It integrates concepts from Single-Cell ProtEomics by Mass Spectrometry (SCoPE-MS) with automated and miniaturized sample preparation, substantially lowering cost and hands-on time. SCoPE2 uses data-driven analytics to optimize instrument parameters for sampling more ion copies per protein, thus supporting quantification with improved count statistics. These advances enabled us to analyze the emergence of cellular heterogeneity as homogeneous monocytes differentiated into macrophage-like cells in the absence of polarizing cytokines. We used SCoPE2 to quantify over 2,000 proteins in 356 single monocytes and macrophages in about 85 hours of instrument time, and the quantified proteins allowed us to discern single cells by cell type. Furthermore, the data uncovered a continuous gradient of proteome states for the macrophage-like cells, suggesting that macrophage heterogeneity may emerge even in the absence of polarizing cytokines. Our methodology lays the foundation for quantitative analysis of protein networks at single-cell resolution.",
              url = "https://doi.org/10.1101/665307",
              dateStamp = "2019-07-09",
              name = "Harrison Specht, Edward Emmott, David H. Perlman, Antonius Koller, Nikolai Slavov",
              lab = "Slavov Lab",
              instrumentModel = "Q-Exactive Orbitrap",
              instrumentManufacturer = "Thermo Scientific",
              softwareName = "MaxQuant",
              softwareVersion = "1.6.2.3",
              switchingCriteria = "After a precursor scan from 450 to 1600 m/z at 70,000 resolving power, the top 5 most intense precursor ions with charges 2 to 4 and above the AGC min threshold of 20,000 were isolated for MS2 analysis via a 0.7 Th isolation window",
              ionSource = "ESI",
              ionSourceDetails = "Electrospray voltage was set to 2,200V, applied at the end of the analytical column. To reduce atmospheric background ions and enhance peptide signal to noise ratio, an Active Background Ion Reduction Device (ABIRD, by ESI Source Solutons, LLC, Woburn MA, USA) was used at the nanospray interface. The temperature of ion transfer tube was 250 degrees Celsius and the S-lens RF level set to 80.",
              analyser = "orbitrap",
              analyserDetails = "Precursor ions were accumulated for at most 300ms. Then they were fragmented via HCD at a and the fragments analyzed at 70,000 resolving power. Dynamic exclusion was used with a duration of 30 seconds with a mass tolerance of 10ppm.",
              collisionEnergy = "33 eV (normalized to m/z 500, z=1)")


####---- Peptide data ----####

## This chunk of code saves the peptide data as provide by the authors, that is 
## already preprocessed to some extend
if(FALSE){ 
  # Load metadata
  pdat <- t(read.csv("../extdata/specht2019/Cells.csv", row.names = 1))
  # Load data
  dat <- read.csv(file = "../extdata/specht2019/Peptides-raw.csv")
  rownames(dat) <- dat[, 2]
  
  # Create the expression data for the MSnSet object
  edat <- as.matrix(dat[, -c(1,2)])
  
  # Create the feature data for the MSnSet object
  fdat <- AnnotatedDataFrame(data = dat[, 1:2])
  
  # Create the MSnSet object
  specht2019_peptide <- new("MSnSet", exprs = edat, 
                            featureData = fdat,
                            phenoData = pdat,
                            experimentData = expdat)
}

## Load the PSM data
dat0 <- read.table("../extdata/specht2019/ev_updated.txt", header = TRUE, sep = "\t")
## Before formatting data, keep only relevant fields 
intensity.coln <- colnames(dat0)[grepl("^Reporter[.]intensity[.]\\d+$", colnames(dat0))]
.keep <- c("Raw.file", "Modified.sequence", "Sequence", "Length", "Charge",
           "Proteins", "Reverse", "Leading.razor.protein", 
           "Potential.contaminant", "Gene.names", "Protein.names", "Mass", 
           "Calibrated.retention.time", "PIF", "dart_PEP",  "dart_qval", 
           "Score", intensity.coln)
dat <- dat0[,.keep]

## Correct for isotopic cross contamination from the TMT-11 channel
if(FALSE){ # This is not needed because MaxQuant did the correction already
  icc <- read.csv("../extdata/specht2019/te269088_lot_correction.csv", 
                  row.names = 1)
  corrected.ri <- t( solve(icc) %*% t(dat[,intensity.coln]) )
  corrected.ri[corrected.ri < 0.1] <- NA
  dat[,intensity.coln] <- corrected.ri
}

## Filter data
# Keep only single cell data and remove experiment FP95 and FP96 (lower quality)
dat <- dat[!grepl("col19|col2[0-4]|blank|QC|FP9[56]", dat$Raw.file), ]
# Remove the reverse hits (from decoy database) and contaminants
dat <- dat[dat$Reverse != "+", ]
dat <- dat[!grepl("^REV", dat$Leading.razor.protein), ]
dat <- dat[dat$Potential.contaminant != "+", ]
dat <- dat[dat$PIF > 0.8 & !is.na(dat$PIF), ]

# Remove spectra with poor identification confidence
# This is done with the output of DART-ID
qprots <- unique(dat[dat$dart_qval < 0.01, "Leading.razor.protein"])
dat <- dat[dat$Leading.razor.protein %in% qprots, ]
dat <- dat[dat$dart_PEP < 0.02, ]

# Remove runs with insufficient peptides 
pep.t <- table(dat$Raw.file)
dat <- dat[! dat$Raw.file %in% names(pep.t[pep.t < 300]),] # 25 runs were removed

## Format data
# Change the table to long format by merging the 11 reporter intensities as 2 
# distinct variables (TMT channel and signal intensity)
dat <- as.data.frame(pivot_longer(data = dat, cols = intensity.coln, 
                                  names_to = "Channel", 
                                  values_to = "Reporter.intensity"))
# Deal with duplicate peptides (measured as differently charged ions)
# In case of duplicate peptides in the same run for the same channel, we keep 
# the  peptide with the lowest PEP
dat <- dat[order(dat$dart_PEP, decreasing = FALSE),]
dat <- dat[!duplicated(dat[,c("Modified.sequence", "Raw.file", "Channel")]), ]
# Create a unique ID for samples 
dat$sample <- paste0(dat$Raw.file, "-", dat$Channel)

## Format the expression data
edat <- pivot_wider(dat, id_cols = "Modified.sequence", names_from = "sample",
                    values_from = "Reporter.intensity", 
                    values_fill = list("Reporter.intensity" = NA))
# 0 intensities are missing values
edat[edat == 0] <- NA
# Add rownames
rown <- edat[,1]
edat <- as.matrix(edat[,-1])
rownames(edat) <- rown[[1]]

# Create the phenotype data 
pdat <- do.call(rbind, lapply(colnames(edat), function(x){
  x <- strsplit(x, "-")[[1]]
  data.frame(run = x[1], channel = x[2])
}))
pdat$channel <- as.numeric(sub("[^\\d]*[.]", "", pdat$channel)) + 1
rownames(pdat) <- colnames(edat)
# Add the sample info
annot <- read.csv("../extdata/specht2019/annotation_fp60-97.csv")
pdat$sampleType <- sapply(1:nrow(pdat), function(i){
  annot[pdat$channel[i], paste0("X", pdat$run[i])]
})
# Add the batch info
batch <- read.csv("../extdata/specht2019/batch_fp60-97.csv", row.names = 1)
pdat <- cbind(pdat, batch[as.character(pdat$run),])

# Create the feature data 
# Keep only variable that are common accross runs
fdat <- dat[, c("Modified.sequence", "Sequence", "Length", "Proteins", "Charge",
              "Gene.names", "Protein.names", "Mass")]
fdat <- fdat[!duplicated(fdat$Modified.sequence), ]
rownames(fdat) <- fdat$Modified.sequence
fdat <- fdat[rownames(edat),]

# Create the MSnSet object
# specht2019_peptide <- new("MSnSet", exprs = edat, 
#                            featureData = fdat,
#                            phenoData = pdat,
#                            experimentData = expdat)
# Create the SingleCellExperiment object
specht2019_peptide <-
  SingleCellExperiment(assay = list(peptide = edat), 
                       rowData = fdat,
                       colData = pdat,
                       metadata = list(experimentData = expdat))
# Save data as Rda file
# Note: saving is assumed to occur in "scpdata/inst/scripts"
save(specht2019v1_peptide, file = file.path("../../data/specht2019_peptide.rda"),
     compress = "xz", compression_level = 9)


####---- Protein data ---####

# Load data
dat <- read.csv(file = "../extdata/specht2019/Proteins-processed.csv", row.names = 1)

# Create the expression data for the MSnSet object
edat <- as.matrix(dat[, -ncol(dat)])

# Create the feature data for the MSnSet object
fdat <- dat[, ncol(dat), drop = FALSE]

# Load sample metadata
pdat <- t(read.csv("../extdata/specht2019/Cells.csv", row.names = 1))
pdat <- as.data.frame(pdat)

# Create the MSnSet object
# specht2019_protein <- new("MSnSet", exprs = edat, 
#                           featureData = AnnotatedDataFrame(data = fdat),
#                           phenoData = AnnotatedDataFrame(data = pdat),
#                           experimentData = expdat)
# Create the SingleCellExperiment object
specht2019_protein <-
  SingleCellExperiment(assay = list(protein = edat), 
                       rowData = fdat,
                       colData = pdat,
                       metadata = list(experimentData = expdat))
# Save data as Rda file
# Note: saving is assumed to occur in "(...)/scpdata/inst/scripts"
save(specht2019v1_protein, file = file.path("../../data/specht2019_protein.rda"),
     compress = "xz", compression_level = 9)










