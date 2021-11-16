
####---- Williams et al. 2020 ---####


## Williams, Sarah M., Andrey V. Liyu, Chia-Feng Tsai, Ronald J. Moore,
## Daniel J. Orton, William B. Chrisler, Matthew J. Gaffrey, et al. 
## 2020. “Automated Coupling of Nanodroplet Sample Preparation with
## Liquid Chromatography-Mass Spectrometry for High-Throughput 
## Single-Cell Proteomics.” Analytical Chemistry 92 (15): 10588–96.

## The data was downloaded from the MASSIVE repo: MSV000085230. 

## This script formats the TMT dataset with AML single cells

library(SingleCellExperiment)
library(scp)
library(tidyverse)
dataDir <-"~/PhD/.localdata/SCP/williams2020/"

####---- Load peptide data ----####

f <- list.files(dataDir, pattern = "Peptides_AML_SingleCell.txt", 
                recursive = TRUE, full.names = TRUE)
peps <- read.delim(f)
## Remove unnecessary columns
peps <- peps[, !grepl("count|^Inten|[yd].\\d*$", colnames(peps))]
## Get RI intensities
pepsIntensity <- peps[, !grepl("corrected", colnames(peps))]
pepsIntensity <- readSingleCellExperiment(pepsIntensity,
                                          ecol = grep("^Repo.*intensity.", colnames(pepsIntensity)),
                                          fnames = "Sequence")
colnames(pepsIntensity) <- gsub("Rep.*intensity.", "", colnames(pepsIntensity))
## Get corrected RI intensities
pepsCorrected <- peps[, !grepl("intensity.\\d", colnames(peps))]
pepsCorrected <- readSingleCellExperiment(pepsCorrected,
                                          ecol = grep("^Repo.*intensity.", colnames(pepsCorrected)),
                                          fnames = "Sequence")
colnames(pepsCorrected) <- gsub("Rep.*corrected.", "", colnames(pepsCorrected))

####---- Load protein data ----####

f <- list.files(dataDir, pattern = "ProteinGroups_AML_SingleCell.txt", 
                recursive = TRUE, full.names = TRUE)
prots <- read.delim(f)
## Remove unnecessary columns
prots <- prots[, !grepl("nique.peptides.|^Peptides.[A-C]|^Identif|^Intensity.|y.count|Sequence.coverage", colnames(prots))]
## Get RI intensities
protsIntensity <- prots[, !grepl("corrected", colnames(prots))]
protsIntensity <- readSingleCellExperiment(protsIntensity,
                                           ecol = grep("^Repo.*intensity.", colnames(protsIntensity)),
                                           fnames = "Protein.IDs")
colnames(protsIntensity) <- gsub("Rep.*intensity.", "", colnames(protsIntensity))
## Get corrected RI intensities
protsCorrected <- prots[, !grepl("intensity.\\d", colnames(prots))]
protsCorrected <- readSingleCellExperiment(protsCorrected,
                                           ecol = grep("^Repo.*intensity.", colnames(protsCorrected)),
                                           fnames = "Protein.IDs")
colnames(protsCorrected) <- gsub("Rep.*corrected.", "", colnames(protsCorrected))

####---- Load sample annotations ----####

f <- list.files(dataDir, pattern = "TMT_Labelling.txt", 
                recursive = TRUE, full.names = TRUE)
annot <- read.table(f, sep = ":", header = FALSE)
colnames(annot) <- c("Channel", "SampleType")
annot$ChannelIndex <- rownames(annot)
## Add sample amounts (from table S2)
annot$Amount <- c("10ng", "0.2ng", "0", rep("1cell", 8))
## Match sample names
runs <- unique(sub("^\\d*[.]?", "", colnames(pepsIntensity)))
annot <- lapply(runs, function(run) {
    rownames(annot) <- paste0(annot$ChannelIndex, ".", run)
    annot
})
annot <- DataFrame(do.call(rbind, annot))
annot$Batch <- sub("^\\d*.", "", rownames(annot))

####---- Create the QFeatures object ----####

el <- ExperimentList(peptides_intensity = pepsIntensity,
                     peptides_corrected = pepsCorrected,
                     proteins_intensity = protsIntensity,
                     proteins_corrected = protsCorrected)
williams2020_tmt <- QFeatures(el, colData = annot)

## Create assay links
williams2020_tmt <- addAssayLink(williams2020_tmt, 
                             from = "peptides_intensity", 
                             to = "proteins_intensity",
                             varFrom = "Leading.razor.protein", 
                             varTo = "Majority.protein.IDs")
williams2020_tmt <- addAssayLink(williams2020_tmt, 
                             from = "peptides_corrected", 
                             to = "proteins_corrected",
                             varFrom = "Leading.razor.protein", 
                             varTo = "Majority.protein.IDs")

### Save data as Rda file
## Note: saving is assumed to occur in "scpdata/inst/scripts"
save(williams2020_tmt, 
     file = file.path("~/PhD/.localdata/scpdata/williams2020_tmt.Rda"),
     compress = "xz", 
     compression_level = 9)
