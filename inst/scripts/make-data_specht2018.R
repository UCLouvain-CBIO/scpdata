
####---- Specht et al. 2018 ---####


# Specht, Harrison, Guillaume Harmange, David H. Perlman, Edward Emmott, Zachary 
# Niziolek, Bogdan Budnik, and Nikolai Slavov. 2018. “Automated Sample 
# Preparation for High-Throughput Single-Cell Proteomics.” bioRxiv. 
# https://doi.org/10.1101/399774.

library(scp)
library(SingleCellExperiment)
library(tidyverse)
setwd("./inst/scripts")
dataDir <- "../extdata/specht2018/"

# The data was downloaded from 
# https://drive.google.com/drive/folders/19YG70I52DH5yETcZagdUjNZWNPs0JXVr
# to scpdata/inst/extdata/specht2018/


####---- Peptide data ----####


# Load the PSM data (MaxQuant evidence file)
paste0(dataDir, "evidence.txt") %>%
  read.table(header = TRUE, sep = "\t") %>%
  rename(Set = Raw.file) ->
  ev

## Create the sample metadata
sets <- unique(ev$Set)
channels <- grep("^R.*y[.][0-9]*$", colnames(ev), value = TRUE)
data.frame(Set = rep(sets, each = length(channels)),
           Channel = rep(channels, length(sets)))
meta


## STOPED HERE: NEED THE SAMPLE ANNOTATION TABLE FOR THIS DATASET TO BE USEFUL!!




dat <- dat0
# Column names holding the TMT intensities
intensity.coln <- colnames(dat)[grepl("^Reporter[.]intensity[.]\\d+$", colnames(dat))]

# Sample data 
# TODO ask for exp design

## Filter data

# Remove the reverse hits (from decoy database) and contaminants
dat <- dat[dat$Reverse != "+", ] # sum(dat0$Reverse == "+")/nrow(dat0) * 100 # 20.02% 
dat <- dat[!grepl("^REV", dat$Leading.razor.protein), ] # sum(grepl("^REV", dat0$Leading.razor.protein))/nrow(dat0) * 100 # 20.02% 
dat <- dat[dat$Potential.contaminant != "+", ] # sum(dat0$Potential.contaminant == "+")/nrow(dat0) * 100 # 5.01% 
dat <- dat[dat$PIF > 0.8 & !is.na(dat$PIF), ] # sum(dat0$PIF <= 0.8 | is.na(dat0$PIF))/nrow(dat0) * 100 # 18.21% 

# Remove spectra with poor identification confidence
# TODO missing FDR ...
dat <- dat[dat$PEP < 0.02, ] # sum(dat0$PEP > 0.02)/nrow(dat0) * 100 # 51.00 % 

# Remove runs with insufficient peptides 
pep.t <- table(dat$Raw.file)
dat <- dat[! dat$Raw.file %in% names(pep.t[pep.t < 300]),]

# Other filters based on carrier ratio
# TODO when experiment design data is available

## Format data

# Change the table to long format by merging the 11 reporter intensities as 2 
# distinct variables (TMT channel and signal intensity)
dat <- pivot_longer(data = dat, cols = intensity.coln, 
                    names_to = "Channel", 
                    values_to = "Reporter.intensity")

# Deal with duplicate peptides (measured as differently charged ions)
# In case of duplicate peptides in the same run for the same channel, we keep the
# peptide with the lowest PEP
dat <- dat[order(dat$PEP, decreasing = FALSE),]
dat <- dat[!duplicated(dat[,c("Modified.sequence", "Raw.file", "Channel")]), ]

# Create a unique ID for samples 
dat$sample <- paste0(dat$Raw.file, "-", dat$Channel)

# Format the expression data
edat <- pivot_wider(dat, id_cols = "Modified.sequence", names_from = "sample",
                    values_from = "Reporter.intensity", 
                    values_fill = list("Reporter.intensity" = NA))
# 0 intensities are missing values
edat[edat == 0] <- NA
# Add rownames
rown <- edat[,1]
edat <- as.matrix(edat[,-1])
rownames(edat) <- rown$Modified.sequence

# Create the phenotype data 
pdat <- do.call(rbind, lapply(colnames(edat), function(x){
  x <- strsplit(x, "-")[[1]]
  data.frame(run = x[1], channel = x[2])
}))
rownames(pdat) <- colnames(edat)
# Create sample metadata description
mpdat <- data.frame(type = c("characters",
                             "characters"),
                    labelDescription = c("the id of the SCoPE-MS run set",
                                         "TMT channel id within the SCoPE-MS set"))
# Combine all phenotype data in an AnnotatedDataFrame
pdat <- AnnotatedDataFrame(data = pdat, varMetadata = mpdat)

# Create the feature data 
fdat <- do.call(rbind, lapply(rownames(edat), function(x){
  .idx <- dat$Modified.sequence == x
  .sub <- dat[.idx, , drop = FALSE]
  # Add common field 
  line <- unique(.sub[, c("Modified.sequence", "Sequence", "Length", "Proteins", 
                          "Gene.names", "Protein.names", "Mass"),])
  # Add peptide count, that is in how many samples was the peptide found
  line$Count <- sum(.idx)
  # Average the remaining fields fields
  line <- data.frame(line, 
                     as.list(apply(.sub[,c("Calibrated.retention.time", "PIF", 
                                           "PEP", "Score")], 2, median, na.rm = TRUE)))
  return(line)
}))
rownames(fdat) <- rownames(edat)
# Make data description for the feature metadata
mfdat <- data.frame(type = c("characters",
                             "characters",
                             "integer",
                             "character",
                             "character",
                             "character",
                             "double",
                             "integer",
                             "double",
                             "double",
                             "double",
                             "double"),
                    labelDescription = c("Sequence representation including the post-translational modifications (abbreviation of the modification in brackets before the modified AA). The sequence is always surrounded by underscore characters ('_').",
                                         "The identified amino acid sequence of the peptide.",
                                         "Number of amino acids of the identified peptide.",
                                         "The identifiers of the proteins a particular peptide is associated with.",
                                         "Name(s) of gene(s) the peptide is associated with.",
                                         "Name(s) of the protein(s) the peptide is associated with.",
                                         "The predicted monoisotopic mass of the identified peptide sequence.",
                                         "The number of samples where the peptide was identified.",
                                         "The recalibrated retention time in minutes in the elution profile of the precursor ion. The value is the median over all samples where the peptide was identified.",
                                         "Short for Parent Ion Fraction; indicates the fraction the target peak makes up of the total intensity in the inclusion window. The value is the median over all samples where the peptide was identified.",
                                         "Posterior Error Probability of the identification. This value essentially operates as a p-value, where smaller is more significant. The value is the median over all samples where the peptide was identified.",
                                         "Andromeda score for the best associated MS/MS spectrum. The value is the median over all samples where the peptide was identified."))
## Combine all feature data in an AnnotatedDataFrame
fdat <- AnnotatedDataFrame(data = fdat, varMetadata = mfdat)

## Create the experimental metadata
expdat <- new("MIAPE",
              title = "Automated sample preparation for high-throughputsingle-cell proteomics",
              abstract = "A major limitation to applying quantitative LC-MS/MS proteomics to small samples, suchas single cells, are the losses incured during sample cleanup. To relieve this limitation, we de-veloped a Minimal ProteOmic sample Preparation (mPOP) method for culture-grown mam-malian cells. mPOP obviates cleanup and thus eliminates cleanup-related losses while expe-diting sample preparation and simplifying its automation.  Bulk SILAC samples processedby mPOP or by conventional urea-based methods indicated that mPOP results in completecell lysis and accurate relative quantification. We integrated mPOP lysis with the Single CellProtEomics by Mass Spectrometry (SCoPE-MS) sample preparation, and benchmarked thequantification  of  such  samples  on  a  Q-exactive  instrument.   The  results  demonstrate  lownoise and high technical reproducibility. Then, we FACS sorted single U-937, HEK-293, andmouse ES cells into 96-well plates and analyzed them by automated mPOP and SCoPE-MS.The quantified proteins enabled separating the single cells by cell-type and cell-division-cyclephase.",
              url = "http://dx.doi.org/10.1101/399774",
              dateStamp = "2018-08-25",
              name = "Harrison Specht, Guillaume Harmange, David H. Perlman, Edward Emmott, Zachary Niziolek, Bogdan Budnik, Nikolai Slavov",
              lab = "Slavov Lab",
              instrumentModel = "Q-Exactive Orbitrap",
              instrumentManufacturer = "Thermo Scientific",
              softwareName = "MaxQuant",
              softwareVersion = "1.6.2.3",
              switchingCriteria = "the top 5 most intense precursors with charges 2 to 4 were selected for HCD fragmentation at resolution 70,000 with a max fill time of 300ms.",
              ionSource = "not specified",
              ionSourceDetails = "not specified",
              analyser = "orbitrap",
              analyserDetails = "a 0.7 Th isolation window was used for MS2 scans.",
              collisionEnergy = "not specified")

## Create the MSnSet object
# specht2018_peptide <- new("MSnSet", exprs = edat, 
#                           featureData = fdat,
#                           phenoData = pdat,
#                           experimentData = expdat)
## Create SingleCellExperiment object
specht2018_peptide <-
  SingleCellExperiment(assay = list(peptide = edat), 
                       rowData = DataFrame(pData(fdat)),
                       colData = DataFrame(pData(pdat)),
                       metadata = list(experimentData = expdat))
## Save data as Rda file
## Note: saving is assumed to occur in "(...)/scpdata/inst/scripts"



