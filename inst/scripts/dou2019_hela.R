
####---- Dou et al. 2019 - HeLa digests ----####


## Dou, Maowei, Geremy Clair, Chia-Feng Tsai, Kerui Xu, William B. Chrisler, 
## Ryan L. Sontag, Rui Zhao, et al. 2019. “High-Throughput Single Cell Proteomics 
## Enabled by Multiplex Isobaric Labeling in a Nanodroplet Sample Preparation 
## Platform.” Analytical Chemistry, September. 
## https://doi.org/10.1021/acs.analchem.9b03349.

## This article contains 3 datasets. The script will focus on the HeLa digests

library(openxlsx)
library(MSnbase)
library(SingleCellExperiment)
setwd("inst/scripts/")


####---- Experimental metadata ----####

## Create the experimental metadata (common to all data sets)
expdat <- new("MIAPE",
              title = "High-Throughput Single Cell Proteomics Enabled by Multiplex Isobaric Labeling in a Nanodroplet Sample Preparation Platform",
              abstract = "Effective extension of mass spectrometry-based proteomics to single cells remains challenging. Herein we combined microfluidic nanodroplet technology with tandem mass tag (TMT) isobaric labeling to significantly improve analysis throughput and proteome coverage for single mammalian cells. Isobaric labeling facilitated multiplex analysis of single cell-sized protein quantities to a depth of ∼1 600 proteins with a median CV of 10.9% and correlation coefficient of 0.98. To demonstrate in-depth high throughput single cell analysis, the platform was applied to measure protein expression in 72 single cells from three murine cell populations (epithelial, immune, and endothelial cells) in <2 days instrument time with over 2 300 proteins identified. Principal component analysis grouped the single cells into three distinct populations based on protein expression with each population characterized by well-known cell-type specific markers. Our platform enables high throughput and unbiased characterization of single cell heterogeneity at the proteome level.",
              url = "https://doi.org/10.1021/acs.analchem.9b03349",
              dateStamp = "2019-09-11",
              name = "Dou, Maowei; Clair, Geremy; Tsai, Chia-Feng; Xu, Kerui; Chrisler, William B; Sontag, Ryan L; Zhao, Rui; Moore, Ronald J; Liu, Tao; Pasa-Tolic, Ljiljana; Smith, Richard D; Shi, Tujin; Adkins, Joshua N; Qian, Wei-Jun; Kelly, Ryan T; Ansong, Charles; Zhu, Ying",
              lab = "Pacific Northwest National Laboratory (PNNL) - Environmental Molecular Sciences Laboratory",
              instrumentModel = "Orbitrap Fusion Lumos Tribrid mass spectrometer",
              instrumentManufacturer = "Thermo Fischer",
              softwareName = "MS-GF+ (identification) and MASIC (quantification)",
              softwareVersion = "not specified",
              switchingCriteria = "Precursor ions with charges of +2 to +7 and intensities >20,000 were selected for MS/MS sequencing.",
              ionSource = "ESI",
              ionSourceDetails = "The peptides were ionized at a spray voltage of 2 kV and the ions were collected into an ion transfer capillary set at 150 °C. The RF lens was set at 30%.",
              analyser = "orbitrap",
              analyserDetails = "The MS1 scan was set at a mass  range  from  375  to  1575,  a  scan  resolution  of  120  k,  an  AGC  target  of  1E6,  and  a maximum  injection  time  of  50  ms.  Precursor  ions  with  charges  of  +2  to  +7  and intensities  >20,000  were  selected  for  MS/MS  sequencing.  Precursor  ions  were  isolated with an m/z window of 0.7 Th and fragmented by high energy dissociation (HCD) set at 35%. Repeat sampling were reduced with an exclusion duration of 60 s and m/z tolerance of ±10 ppm. MS/MS scan was carried out in the Orbitrap with an AGC target of 1E5. The MS/MS scan resolutions and maximum injection times were set as 60 k and 246 ms.",
              collisionEnergy = "HCD of 35%")


####---- Peptide data ----####

## Load and combine the identification and quantification data
expNames <- c("Hela_run_1", "Hela_run_2")
## For every experiment (= MS run)
dats <- lapply(expNames, function(exp){
  ## Identification data
  idFile <- list.files("../../extdata/dou2019/mzid/", pattern = exp, 
                  full.names = TRUE)
  mzidFile <- openIDfile(idFile)
  mzid <- psms(mzidFile)
  mzid <- mzid[!duplicated(mzid), ]
  ## Remove PSMs that match decoy database
  mzid <- mzid[!mzid$isDecoy, ]
  ## Remove contaminants
  mzid <- mzid[!grepl("Contaminant", mzid$DatabaseAccess), ]
  
  ## Quantification data
  quantFile <- list.files("../../extdata/dou2019/quant/", pattern = exp, 
                  full.names = TRUE)
  quant <- read.table(quantFile, header = TRUE, sep = "\t")
  ## Replace 0's with NA's
  quant[quant == 0] <- NA
  
  ## Combined data 
  matchInds <- match(x = mzid$scan.number.s., 
                     table = quant$ScanNumber, 
                     nomatch = NA)
  dat <- cbind(mzid, quant[matchInds, ])
  
  ## Create the MSnSet object
  ecols <- grepl("^Ion_.*\\d$", colnames(dat))
  ## TODO should scpdata also provide PSM data -> dou2019_hela_psm ? 
  dat <- psmDat <- MSnSet(exprs = as.matrix(dat[, ecols]), 
                          fData = dat[, !ecols])
  
  ## Perform isotope correction
  impurities <- makeImpuritiesMatrix(10, edit = FALSE) # TMT 10 was used
  dat <- purityCorrect(dat, impurities)
  
  ## Group the observation by peptide (modified and charged)
  dat <- combineFeatures(object = dat, 
                         groupBy = fData(dat)$sequence, 
                         method = "sum", ## As done at protein level (see doc)
                         cv = FALSE, 
                         na.rm = TRUE)
  
  ## Filter only relevant feature data vars
  keep <- c(peptide = "sequence", 
            protein = "DatabaseAccess", 
            proteinLength = "DBseqLength", 
            proteinDescr = "DatabaseDescription",
            start = "start", 
            end = "end")
  fData(dat) <- fData(dat)[, keep]
  fvarLabels(dat) <- names(keep) ## Replace labels with more intuitive naming
  
  ## Extract the phenotype data
  pData(dat) <- data.frame(Ion = gsub("Ion_", "", colnames(dat)),
                           ExperimentLabel = exp,
                           Run = as.factor(sub("Hela_run_", "", exp)), 
                           row.names = colnames(dat))
  dat <- updateSampleNames(dat, label = exp, sep = "_")
  
  ## Return data
  cat("Formated:", exp, "\n")
  return(dat)
})

## Merge the different runs
dat <- do.call(combine, dats)
## Note warnings are thrown because the levels in pData and fData don't match 
## between runs. 

## Add experimental data
experimentData(dat) <- expdat

## Create SingleCellExperiment object
dou2019_hela_peptide <- 
  SingleCellExperiment(assay = SimpleList(peptide = exprs(dat)),
                       rowData = data.frame(fData(dat)), 
                       colData = data.frame(pData(dat)),
                       metadata = list(experimentData = expdat))
## Save data as Rda file
## Note: saving is assumed to occur in "scpdata/inst/scripts"
save(dou2019_hela_peptide, compress = "xz", compression_level = 9,
     file = file.path("../../data/dou2019_hela_peptide.rda"))


####---- Protein data ----####


## Load the data
## Data was downloaded from https://doi.org/10.1021/acs.analchem.9b03349.
dataFile <- "../../extdata/dou2019/protein/ac9b03349_si_003.xlsx"
dat_xlsx <- loadWorkbook(dataFile)
dat <- read.xlsx(dat_xlsx, sheet = 6, colNames = FALSE)

## Extract the feature data
fd <- dat[, 1:3]
colnames(fd) <- fd[3, ]
fd <- fd[-(1:3), ]
rownames(fd) <- fd$Protein
fd[, 2:3] <- sapply(fd[, 2:3], as.numeric)
fd <- DataFrame(fd)

## Extract the phenotype data
run <- sub("run", "", unlist(dat[1, -(1:3)]))
sampleType <- paste0(unlist(dat[2, -(1:3)]), c(1:7, 1:2, ""))
TMT <- sapply(dat[3, -(1:3)], function(x) strsplit(x, "_")[[1]][2])
datasetID <- sapply(dat[3, -(1:3)], function(x) strsplit(x, "_")[[1]][4])
pd <- DataFrame(sampleType, run, TMT, datasetID)
rownames(pd) <- paste0(sampleType, "_", run)

## Extract the expression data
ed <- apply(dat[-(1:3), -(1:3)], 2, as.numeric)
rownames(ed) <- rownames(fd)
colnames(ed) <- rownames(pd)

## Create MSnSet object
# dou2019_hela_protein <- new("MSnSet", exprs = ed,
#                             featureData = AnnotatedDataFrame(fd),
#                             phenoData = AnnotatedDataFrame(pd),
#                             experimentData = expdat)

## Create SingleCellExperiment object
dou2019_hela_protein <- 
  SingleCellExperiment(assay = SimpleList(protein = ed),
                       rowData = fd, 
                       colData = pd,
                       metadata = list(experimentData = expdat))
## Save data as Rda file
## Note: saving is assumed to occur in "scpdata/inst/scripts"
save(dou2019_hela_protein, compress = "xz", compression_level = 9,
     file = file.path("../../data/dou2019_hela_protein.rda"))





