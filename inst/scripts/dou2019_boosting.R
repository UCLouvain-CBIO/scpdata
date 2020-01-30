
####---- Dou et al. 2019 - testing boosting ratios ----####


## Dou, Maowei, Geremy Clair, Chia-Feng Tsai, Kerui Xu, William B. Chrisler, 
## Ryan L. Sontag, Rui Zhao, et al. 2019. “High-Throughput Single Cell Proteomics 
## Enabled by Multiplex Isobaric Labeling in a Nanodroplet Sample Preparation 
## Platform.” Analytical Chemistry, September. 
## https://doi.org/10.1021/acs.analchem.9b03349.

## This article contains 3 datasets. The script will focus on the testing of
## boosting ratios 

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
              softwareName = "MS-GF+",
              softwareVersion = "not specified",
              switchingCriteria = "Precursor ions with charges of +2 to +7 and intensities >20,000 were selected for MS/MS sequencing.",
              ionSource = "ESI",
              ionSourceDetails = "The peptides were ionized at a spray voltage of 2 kV and the ions were collected into an ion transfer capillary set at 150 °C. The RF lens was set at 30%.",
              analyser = "orbitrap",
              analyserDetails = "The MS1 scan was set at a mass  range  from  375  to  1575,  a  scan  resolution  of  120  k,  an  AGC  target  of  1E6,  and  a maximum  injection  time  of  50  ms.  Precursor  ions  with  charges  of  +2  to  +7  and intensities  >20,000  were  selected  for  MS/MS  sequencing.  Precursor  ions  were  isolated with an m/z window of 0.7 Th and fragmented by high energy dissociation (HCD) set at 35%. Repeat sampling were reduced with an exclusion duration of 60 s and m/z tolerance of ±10 ppm. MS/MS scan was carried out in the Orbitrap with an AGC target of 1E5. The MS/MS scan resolutions and maximum injection times were set as 60 k and 246 ms.",
              collisionEnergy = "HCD of 35%")

## Sample annotation provided in Table S2 in the Supplementary information
annot <- data.frame(Ion126 = c("Raw", "SVEC", "C10", "SVEC", "C10", "SVEC"),
                    Ion127N = c("Raw", "SVEC", "C10", "SVEC", "C10", "SVEC"),
                    Ion127C = c("Raw", "SVEC", "Raw", "SVEC", "C10", "SVEC"),
                    Ion128N = c("SVEC", "C10", "Raw", "C10", "Raw", "C10"),
                    Ion128C = c("SVEC", "C10", "Raw", "C10", "Raw", "C10"),
                    Ion129N = c("SVEC", "C10", "SVEC", "C10", "Raw", "C10"),
                    Ion129C = c("C10", "Raw", "SVEC", "Raw", "SVEC", "Raw"),
                    Ion130N = c("C10", "Raw", "empty", "empty", "empty", "empty"),
                    Ion130C = c("C10", "Raw", "C10", "Raw", "SVEC", "Raw"),
                    Ion131N = c("empty", "empty", "boost", "boost", "boost", "boost"), 
                    row.names = paste0("Boosting_", rep(c(0, 5, 50), each = 2), 
                                       "ng_run_", rep(1:2, 3)),
                    stringsAsFactors = FALSE)


####---- Peptide data ----####

## Load and combine the identification and quantification data
expNames <- list.files("../../extdata/dou2019/mzid/")
expNames <- grep("^Boosting", expNames, value = TRUE)
expNames <- sub("(Boosting.*)_msgf.*$", "\\1", expNames)
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
  
  ## Combined data 
  matchInds <- match(x = mzid$scan.number.s., 
                     table = quant$ScanNumber, 
                     nomatch = NA)
  dat <- cbind(mzid, quant[matchInds, ])
  
  ## Create the MSnSet object
  ecols <- grepl("^Ion_.*\\d$", colnames(dat))
  ## TODO should scpdata also provide PSM data -> dou2019_mouse_psm ? 
  dat <- psmDat <- MSnSet(exprs = as.matrix(dat[, ecols]), 
                          fData = dat[, !ecols])
  
  ## Perform isotope correction
  impurities <- makeImpuritiesMatrix(10, edit = FALSE) # TMT 10 was used
  dat <- purityCorrect(dat, as.matrix(impurities))
  
  
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
            startPos = "start", 
            endPos = "end")
  fData(dat) <- fData(dat)[, keep]
  fvarLabels(dat) <- names(keep) ## Replace labels with more intuitive naming
  
  ## Extract the phenotype data
  pData(dat) <- data.frame(Ion = sub("Ion_", "", colnames(dat)),
                           ExperimentLabel = exp,
                           SampleType = unlist(annot[exp, ]),
                           Run = sub("^.*run_(\\d)$", "\\1", exp),
                           Boosting = sub("^.*ing_(\\d*)ng.*$", "\\1", exp),
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
## Replace NA's by 0's otherwise some SingleCellExperiment functions break
## TODO this seems very dirty to me...
exprs(dat)[is.na(exprs(dat))] <- 0 

## Create SingleCellExperiment object
dou2019_boosting_peptide <- 
  SingleCellExperiment(assay = SimpleList(peptide = exprs(dat)),
                       rowData = data.frame(fData(dat)), 
                       colData = data.frame(pData(dat)),
                       metadata = list(experimentData = expdat))
## Save data as Rda file
## Note: saving is assumed to occur in "scpdata/inst/scripts"
save(dou2019_boosting_peptide, compress = "xz", compression_level = 9,
     file = file.path("../../data/dou2019_boosting_peptide.rda"))


####---- Protein data ----####

## Load the data
## Data was downloaded from https://doi.org/10.1021/acs.analchem.9b03349.
dataFile <- "../../extdata/dou2019/protein/ac9b03349_si_004.xlsx"
dat_xlsx <- loadWorkbook(dataFile)
dat_nb <- read.xlsx(dat_xlsx, sheet = 2, colNames = FALSE)
dat_5b <- read.xlsx(dat_xlsx, sheet = 4, colNames = FALSE)
dat_50b <- read.xlsx(dat_xlsx, sheet = 6, colNames = FALSE)

## Extract the phenotype data
## Internal function
extractPhenoData <- function(dat, boosting){
  run <- rep(paste0("run", 1:2), each = 10)
  sampleType <- sapply(unlist(dat[1, -(1:3)]), function(x) strsplit(x, "_")[[1]][1])
  TMT <- sapply(dat[2, -(1:3)], function(x) strsplit(x, "_")[[1]][2])
  datasetID <- sapply(dat[2, -(1:3)], function(x) strsplit(x, "_")[[1]][4])
  pd <- data.frame(sampleType, run, TMT, datasetID, boosting)
  rownames(pd) <- make.unique(paste0(sampleType, "_", run, "_", boosting), sep = "_")
  return(pd)
}
## No boosting
pd_nb <- extractPhenoData(dat_nb, "0ng")
## 5ng boosting 
pd_5b <- extractPhenoData(dat_5b, "5ng")
## 50ng boosting 
pd_50b <- extractPhenoData(dat_50b, "50ng")
## Combine phenotype data
pd <- rbind(pd_nb, pd_5b, pd_50b)

## Extract the feature data
extractFeatData <- function(dat){
  fd <- dat[, 1:3]
  colnames(fd) <- fd[2, ]
  fd <- fd[-(1:2), ]
  rownames(fd) <- fd$Protein
  fd[, 2:3] <- sapply(fd[, 2:3], as.numeric)
  return(fd)
}
## No boosting
fd_nb <- extractFeatData(dat_nb)
## 5ng boosting 
fd_5b <- extractFeatData(dat_5b)
## 50ng boosting
fd_50b <- extractFeatData(dat_50b)
## Combine feature data
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

## Extract the expression data
ed_nb <- as.matrix(dat_nb[-(1:2), -(1:3)])
ed_5b <- as.matrix(dat_5b[-(1:2), -(1:3)])
ed_50b <- as.matrix(dat_50b[-(1:2), -(1:3)])
## Combine expression data 
ed <- matrix(NA, nrow = nrow(fd), ncol = nrow(pd),
             dimnames = list(rownames(fd), rownames(pd)))
ed[match(fd_nb$Protein, prots, nomatch = 0), 1:20] <- as.numeric(ed_nb)
ed[match(fd_5b$Protein, prots, nomatch = 0), 21:40] <- as.numeric(ed_5b)
ed[match(fd_50b$Protein, prots, nomatch = 0), 41:60] <- as.numeric(ed_50b)

# ## Create MSnSet object
# dou2019_boosting_protein <- new("MSnSet", exprs = ed,
#                                 featureData = AnnotatedDataFrame(fd),
#                                 phenoData = AnnotatedDataFrame(pd),
#                                 experimentData = expdat)
## Create SingleCellExperiment object
dou2019_boosting_protein <- 
  SingleCellExperiment(assay = SimpleList(protein = ed),
                       rowData = fd, 
                       colData = pd,
                       metadata = list(experimentData = expdat))
## Save data as Rda file
## Note: saving is assumed to occur in "scpdata/inst/scripts"
save(dou2019_boosting_protein, file = file.path("../../data/dou2019_boosting_protein.rda"),
     compress = "xz", compression_level = 9)
