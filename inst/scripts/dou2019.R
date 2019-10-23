
####---- Dou et al. 2019 ----####


# Dou, Maowei, Geremy Clair, Chia-Feng Tsai, Kerui Xu, William B. Chrisler, 
# Ryan L. Sontag, Rui Zhao, et al. 2019. “High-Throughput Single Cell Proteomics 
# Enabled by Multiplex Isobaric Labeling in a Nanodroplet Sample Preparation 
# Platform.” Analytical Chemistry, September. 
# https://doi.org/10.1021/acs.analchem.9b03349.

# This article contains 3 datasets:
# dou2019_1: Supplementary data set 1, raw data for HeLa digest (XLSX)
# dou2019_2r: Supplementary data set 2, raw data for testing boosting ratios (XLSX)
# dou2019_3: Supplementary data set 3, raw data for isobaric labelling-based single cell quantification and bulk-scale label free quantification (XLSX)

library(openxlsx)
library(MSnbase)

### dou2019_1

# Download the data
dataFile <- "https://pubs.acs.org/doi/suppl/10.1021/acs.analchem.9b03349/suppl_file/ac9b03349_si_003.xlsx"
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
welltype <- paste0(unlist(dat[2, -(1:3)]), c(1:7, 1:2, ""))
TMT_ion <- sapply(dat[3, -(1:3)], function(x) strsplit(x, "_")[[1]][2])
dataset_id <- sapply(dat[3, -(1:3)], function(x) strsplit(x, "_")[[1]][4])
pd <- data.frame(welltype, run, TMT_ion, dataset_id)
rownames(pd) <- paste0(welltype, "_", run)

# Extract the expression data
ed <- apply(dat[-(1:3), -(1:3)], 2, as.numeric)
rownames(ed) <- rownames(fd)
colnames(ed) <- rownames(pd$celltype)

# Create the MIAPE
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

# Create MSnSet object
x <- MSnSet(exprs = ed, fData = fd, pData = pd)
experimentData(x) <- expdat

# Save data as Rda file
stopifnot(validObject(x))
assign("dou2019_1", x)
save(dou2019_1, file = file.path("../../data/dou2019_1.rda"),
     compress = "xz", compression_level = 9)
# Note: saving is assumed to occur in "(...)/scpdata/inst/scripts"


