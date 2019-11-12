
####---- Dou et al. 2019 ----####


# Dou, Maowei, Geremy Clair, Chia-Feng Tsai, Kerui Xu, William B. Chrisler, 
# Ryan L. Sontag, Rui Zhao, et al. 2019. “High-Throughput Single Cell Proteomics 
# Enabled by Multiplex Isobaric Labeling in a Nanodroplet Sample Preparation 
# Platform.” Analytical Chemistry, September. 
# https://doi.org/10.1021/acs.analchem.9b03349.

# This article contains 3 datasets that was downloaded from 
# https://doi.org/10.1021/acs.analchem.9b03349 to scpdata/inst/extdata/dou2019
# dou2019_1: Supplementary data set 1, raw data for HeLa digest (XLSX)
# dou2019_2: Supplementary data set 2, raw data for testing boosting ratios (XLSX)
# dou2019_3: Supplementary data set 3, raw data for isobaric labelling-based single cell quantification and bulk-scale label free quantification (XLSX)

library(openxlsx)
library(MSnbase)
# setwd("scpdata/inst/scripts/")


### Experimental data

# Create the experimental metadata (common to all data sets)
expdat <- new("MIAPE",
              title = "High-Throughput Single Cell Proteomics Enabled by Multiplex Isobaric Labeling in a Nanodroplet Sample Preparation Platform",
              abstract = "Effective extension of mass spectrometry-based proteomics to single cells remains challenging. Herein we combined microfluidic nanodroplet technology with tandem mass tag (TMT) isobaric labeling to significantly improve analysis throughput and proteome coverage for single mammalian cells. Isobaric labeling facilitated multiplex analysis of single cell-sized protein quantities to a depth of ∼1 600 proteins with a median CV of 10.9% and correlation coefficient of 0.98. To demonstrate in-depth high throughput single cell analysis, the platform was applied to measure protein expression in 72 single cells from three murine cell populations (epithelial, immune, and endothelial cells) in <2 days instrument time with over 2 300 proteins identified. Principal component analysis grouped the single cells into three distinct populations based on protein expression with each population characterized by well-known cell-type specific markers. Our platform enables high throughput and unbiased characterization of single cell heterogeneity at the proteome level.",
              url = "https://doi.org/10.1021/acs.analchem.9b03349",
              dateStamp = "2019-09-11",
              name = "Dou, Maowei; Clair, Geremy; Tsai, Chia-Feng; Xu, Kerui; Chrisler, William B; Sontag, Ryan L; Zhao, Rui; Moore, Ronald J; Liu, Tao; Pasa-Tolic, Ljiljana; Smith, Richard D; Shi, Tujin; Adkins, Joshua N; Qian, Wei-Jun; Kelly, Ryan T; Ansong, Charles; Zhu, Ying",
              lab = "Pacific Northwest National Laboratory (PNNL) - Environmental Molecular Sciences Laboratory",
              instrumentModel = "Orbitrap Fusion Lumos Tribrid mass spectrometer",
              instrumentManufacturer = "Thermo Fischer",
              softwareName = "MS-GF+",
              softwareVersion = "not specified",
              switchingCriteria = "Precursor ions with charges of +2 to +7 and intensities >20,000 were selected for MS/MS sequencing.",
              ionSource = "ESI",
              ionSourceDetails = "The peptides were ionized at a spray voltage of 2 kV and the ions were collected into an ion transfer capillary set at 150 °C. The RF lens was set at 30%.",
              analyser = "orbitrap",
              analyserDetails = "The MS1 scan was set at a mass  range  from  375  to  1575,  a  scan  resolution  of  120  k,  an  AGC  target  of  1E6,  and  a maximum  injection  time  of  50  ms.  Precursor  ions  with  charges  of  +2  to  +7  and intensities  >20,000  were  selected  for  MS/MS  sequencing.  Precursor  ions  were  isolated with an m/z window of 0.7 Th and fragmented by high energy dissociation (HCD) set at 35%. Repeat sampling were reduced with an exclusion duration of 60 s and m/z tolerance of ±10 ppm. MS/MS scan was carried out in the Orbitrap with an AGC target of 1E5. The MS/MS scan resolutions and maximum injection times were set as 60 k and 246 ms.",
              collisionEnergy = "HCD of 35%")


####---- dou2019_1 ----####


# Download the data
dataFile <- "../extdata/dou2019/ac9b03349_si_003.xlsx"
dat_xlsx <- loadWorkbook(dataFile)
dat <- read.xlsx(dat_xlsx, sheet = 6, colNames = FALSE)

# Extract the feature data
fd <- dat[, 1:3]
colnames(fd) <- fd[3, ]
fd <- fd[-(1:3), ]
rownames(fd) <- fd$Protein
fd[, 2:3] <- sapply(fd[, 2:3], as.numeric)

# Extract the phenotype data
run <- unlist(dat[1, -(1:3)])
sample_type <- paste0(unlist(dat[2, -(1:3)]), c(1:7, 1:2, ""))
TMT_ion <- sapply(dat[3, -(1:3)], function(x) strsplit(x, "_")[[1]][2])
dataset_id <- sapply(dat[3, -(1:3)], function(x) strsplit(x, "_")[[1]][4])
pd <- data.frame(sample_type, run, TMT_ion, dataset_id)
rownames(pd) <- paste0(sample_type, "_", run)

# Extract the expression data
ed <- apply(dat[-(1:3), -(1:3)], 2, as.numeric)
rownames(ed) <- rownames(fd)
colnames(ed) <- rownames(pd$sample_type)

# Create MSnSet object
dou2019_1_protein <- new("MSnSet", exprs = ed, 
                         featureData = AnnotatedDataFrame(fd),
                         phenoData = AnnotatedDataFrame(pd),
                         experimentData = expdat)

# Save data as Rda file
# Note: saving is assumed to occur in "scpdata/inst/scripts"
save(dou2019_1_protein, file = file.path("../../data/dou2019_1_protein.rda"),
     compress = "xz", compression_level = 9)


####---- dou2019_2 ----####


# Clear environment
rm(list = ls())  # but reassign expdat at beginning of script

# Download the data
dataFile <- "../extdata/dou2019/ac9b03349_si_004.xlsx"
dat_xlsx <- loadWorkbook(dataFile)
dat_nb <- read.xlsx(dat_xlsx, sheet = 2, colNames = FALSE)
dat_5b <- read.xlsx(dat_xlsx, sheet = 4, colNames = FALSE)
dat_50b <- read.xlsx(dat_xlsx, sheet = 6, colNames = FALSE)

# Extract the feature data
extractFeatData <- function(dat){
  fd <- dat[, 1:3]
  colnames(fd) <- fd[2, ]
  fd <- fd[-(1:2), ]
  rownames(fd) <- fd$Protein
  fd[, 2:3] <- sapply(fd[, 2:3], as.numeric)
  return(fd)
}
# No boosting
fd_nb <- extractFeatData(dat_nb)
# 5ng boosting 
fd_5b <- extractFeatData(dat_5b)
# 50ng boosting
fd_50b <- extractFeatData(dat_50b)

# Extract the phenotype data
# Internal function
extractPhenoData <- function(dat, boosting){
  run <- rep(paste0("run", 1:2), each = 10)
  sample_type <- sapply(unlist(dat[1, -(1:3)]), function(x) strsplit(x, "_")[[1]][1])
  TMT_ion <- sapply(dat[2, -(1:3)], function(x) strsplit(x, "_")[[1]][2])
  dataset_id <- sapply(dat[2, -(1:3)], function(x) strsplit(x, "_")[[1]][4])
  pd <- data.frame(sample_type, run, TMT_ion, dataset_id, boosting)
  rownames(pd) <- make.unique(paste0(sample_type, "_", run, "_", boosting), sep = "_")
  return(pd)
}
# No boosting
pd_nb <- extractPhenoData(dat_nb, "0ng")
# 5ng boosting 
pd_5b <- extractPhenoData(dat_5b, "5ng")
# 50ng boosting 
pd_50b <- extractPhenoData(dat_50b, "50ng")

# Extract the expression data
ed_nb <- as.matrix(dat_nb[-(1:2), -(1:3)])
ed_5b <- as.matrix(dat_5b[-(1:2), -(1:3)])
ed_50b <- as.matrix(dat_50b[-(1:2), -(1:3)])

# Combine feature data
pd <- rbind(pd_nb, pd_5b, pd_50b)

# Combine phenotype data
prots <- unique(c(fd_nb$Protein, fd_5b$Protein, fd_50b$Protein))
colnams <- as.vector(outer(c("Peptide_count", "Spectra_count"), 
                           c("0ng", "5ng", "50ng"), 
                           paste, sep ="_"))
fd <- matrix(NA, nrow = length(prots), ncol = 6, 
             dimnames = list(prots, colnams))
fd[match(fd_nb$Protein, prots, nomatch = 0), 1:2] <- as.matrix(fd_nb[, 2:3])
fd[match(fd_5b$Protein, prots, nomatch = 0), 3:4] <- as.matrix(fd_5b[, 2:3])
fd[match(fd_50b$Protein, prots, nomatch = 0), 5:6] <- as.matrix(fd_50b[, 2:3])
fd <- data.frame(Protein = prots, fd)

# Combine expression data 
ed <- matrix(NA, nrow = nrow(fd), ncol = nrow(pd),
             dimnames = list(rownames(fd), rownames(pd)))
ed[match(fd_nb$Protein, prots, nomatch = 0), 1:20] <- as.numeric(ed_nb)
ed[match(fd_5b$Protein, prots, nomatch = 0), 21:40] <- as.numeric(ed_5b)
ed[match(fd_50b$Protein, prots, nomatch = 0), 41:60] <- as.numeric(ed_50b)

# Create MSnSet object
dou2019_2_protein <- new("MSnSet", exprs = ed, 
                         featureData = AnnotatedDataFrame(fd),
                         phenoData = AnnotatedDataFrame(pd),
                         experimentData = expdat)

# Save data as Rda file
# Note: saving is assumed to occur in "scpdata/inst/scripts"
save(dou2019_2_protein, file = file.path("../../data/dou2019_2_protein.rda"),
     compress = "xz", compression_level = 9)


####---- dou2019_3 ----####


# Clear environment
rm(list = ls()) # but reassign expdat at beginning of script

# Download the data
dataFile <- "../extdata/dou2019/ac9b03349_si_005.xlsx"
dat_xlsx <- loadWorkbook(dataFile)
dat <- read.xlsx(dat_xlsx, sheet = 2, colNames = FALSE)

# Extract the feature data
fd <- dat[, 1:3]
colnames(fd) <- fd[2, ]
fd <- fd[-(1:2), ]
rownames(fd) <- fd$Protein
fd[, 2:3] <- sapply(fd[, 2:3], as.numeric)

# Extract the phenotype data
sample_type <- unlist(dat[1, -(1:3)])
TMT_ion <- sapply(dat[2, -(1:3)], function(x) strsplit(x, "_")[[1]][2])
dataset_id <- sapply(dat[2, -(1:3)], function(x) strsplit(x, "_")[[1]][4])
pd <- data.frame(sample_type, dataset_id, TMT_ion)
rownames(pd) <- paste0(sample_type, "_", dataset_id, "_", TMT_ion)

# Extract the expression data
ed <- apply(dat[-(1:2), -(1:3)], 2, as.numeric)
rownames(ed) <- rownames(fd)
colnames(ed) <- rownames(pd$sample_type)

# Create MSnSet object
dou2019_3_protein <- new("MSnSet", exprs = ed, 
                         featureData = AnnotatedDataFrame(fd),
                         phenoData = AnnotatedDataFrame(pd),
                         experimentData = expdat)

# Save data as Rda file
# Note: saving is assumed to occur in "scpdata/inst/scripts"
save(dou2019_3_protein, file = file.path("../../data/dou2019_3_protein.rda"),
     compress = "xz", compression_level = 9)





