
####---- Woo et al. 2022 ---####


## Woo, Jongmin, Geremy C. Clair, Sarah M. Williams, Song Feng, 
## Chia-Feng Tsai, Ronald J. Moore, William B. Chrisler, et al. 2022.
## “Three-Dimensional Feature Matching Improves Coverage for 
## Single-Cell Proteomics Based on Ion Mobility Filtering.” Cell 
## Systems 0 (0). https://doi.org/10.1016/j.cels.2022.02.003.

## The data was downloaded from the MASSIVE repo: MSV000085937 

## This script formats the macrophage actiation dataset

library(SingleCellExperiment)
library(scp)
library(tidyverse)
dataDir <-"~/PhD/.localdata/SCP/woo2022/RAW_LPS_SingleCellProteomics/"

####---- Load peptide data ----####

f <- list.files(dataDir, pattern = "^pep", 
                recursive = TRUE, full.names = TRUE)
peps <- read.delim(f)

## Remove unnecessary columns
peps <- peps[, !grepl("^Identif|^Experiment", colnames(peps))]
## Get "raw" intensities
pepsIntensity <- peps[, !grepl("^LFQ|Intens.*Lib$", colnames(peps))]
pepsIntensity <- readSingleCellExperiment(
    pepsIntensity,
    ecol = grep("^Intensity.", colnames(pepsIntensity)),
    fnames = "Sequence"
)
colnames(pepsIntensity) <- gsub("Intensity.", "", colnames(pepsIntensity))
## Get LFQ normalized intensities
pepsLfq <- peps[, !grepl("^Intensity.\\d", colnames(peps))]
pepsLfq <- readSingleCellExperiment(
    pepsLfq,
    ecol = grep("^LFQ", colnames(pepsLfq)),
    fnames = "Sequence"
)
colnames(pepsLfq) <- gsub("LFQ.intensity.", "", colnames(pepsLfq))

####---- Load protein data ----####

f <- list.files(dataDir, pattern = "^prot", 
                recursive = TRUE, full.names = TRUE)
prots <- read.delim(f)
## Remove unnecessary columns
prots <- prots[, !grepl("^Identif|^MS.MS.count|^Peptides[.]|nique.peptides[.]|^Sequence.coverage", colnames(prots))]
## Get intensities
protsIntensity <- prots[, !grepl("^iBAQ|^LFQ|Intens.*Lib$", colnames(prots))]
protsIntensity <- readSingleCellExperiment(
    protsIntensity,
    ecol = grep("^Intensity.", colnames(protsIntensity)),
    fnames = "Protein.IDs"
)
colnames(protsIntensity) <- gsub("Intensity.", "", colnames(protsIntensity))
## Get iBAQ normalized intensities
protsIbaq <- prots[, !grepl("^Intensity.|^LFQ|iBAQ.*Lib$", colnames(prots))]
protsIbaq <- readSingleCellExperiment(
    protsIbaq,
    ecol = grep("^iBAQ.", colnames(protsIbaq)),
    fnames = "Protein.IDs"
)
colnames(protsIbaq) <- gsub("iBAQ.", "", colnames(protsIbaq))
## Get LFQ normalized intensities
protsLfq <- prots[, !grepl("^Intensity.|^iBAQ", colnames(prots))]
protsLfq <- readSingleCellExperiment(
    protsLfq,
    ecol = grep("^LFQ.", colnames(protsLfq)),
    fnames = "Protein.IDs"
)
colnames(protsLfq) <- gsub("LFQ.intensity.", "", colnames(protsLfq))

####---- Get sample annotations ----####

f <- list.files(dataDir, pattern = "^summary", 
                recursive = TRUE, full.names = TRUE)
summary <- read.delim(f)

## Guess sample annotations from the spectra file name
annot <- DataFrame(row.names = colnames(protsIntensity),
                   FileName = summary$Raw.file[match(colnames(protsIntensity), summary$Experiment)])
annot$RunIndex <- sub("^(\\d*)_.*$", "\\1", annot$FileName)
annot$Chip <- sub("^.*Chip(\\d*?)_.*$", "\\1", annot$FileName)
annot$Well <- sub("^.*(?:Chip.|Library_.*?)_(.*?)_.*$", "\\1", annot$FileName)
annot$WellRow <- sub("^.", "", annot$Well)
annot$WellCol <- sub("\\d{2}", "", annot$Well)
annot$Treatment <- sub("^(.*?)_.*$", "\\1", rownames(annot))

####---- Create the QFeatures object ----####

el <- ExperimentList(peptides_intensity = pepsIntensity,
                     peptides_LFQ = pepsLfq,
                     proteins_intensity = protsIntensity,
                     proteins_iBAQ = protsIbaq,
                     proteins_LFQ = protsLfq)
woo2022_macrophage <- QFeatures(el, colData = annot)

## Create assay links
woo2022_macrophage <- addAssayLinkOneToOne(woo2022_macrophage, 
                                           from = "peptides_intensity", 
                                           to = "peptides_LFQ")
woo2022_macrophage <- addAssayLink(woo2022_macrophage, 
                             from = "peptides_intensity", 
                             to = "proteins_intensity",
                             varFrom = "Leading.razor.protein", 
                             varTo = "Majority.protein.IDs")
woo2022_macrophage <- addAssayLink(woo2022_macrophage, 
                             from = "peptides_intensity", 
                             to = "proteins_iBAQ",
                             varFrom = "Leading.razor.protein", 
                             varTo = "Majority.protein.IDs")
woo2022_macrophage <- addAssayLink(woo2022_macrophage, 
                             from = "peptides_intensity", 
                             to = "proteins_LFQ",
                             varFrom = "Leading.razor.protein", 
                             varTo = "Majority.protein.IDs")
plot(woo2022_macrophage)

### Save data as Rda file
save(woo2022_macrophage,
     file = file.path("~/PhD/.localdata/scpdata/woo2022_macrophage.Rda"),
     compress = "xz", 
     compression_level = 9)

