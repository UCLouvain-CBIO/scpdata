
####---- Williams et al. 2020 ---####


## Williams, Sarah M., Andrey V. Liyu, Chia-Feng Tsai, Ronald J. Moore,
## Daniel J. Orton, William B. Chrisler, Matthew J. Gaffrey, et al. 
## 2020. “Automated Coupling of Nanodroplet Sample Preparation with
## Liquid Chromatography-Mass Spectrometry for High-Throughput 
## Single-Cell Proteomics.” Analytical Chemistry 92 (15): 10588–96.

## The data was downloaded from the MASSIVE repo: MSV000085230. 

## This script formats the LFQ dataset with MCF10A single cells

library(SingleCellExperiment)
library(scp)
library(tidyverse)
dataDir <-"~/PhD/.localdata/SCP/williams2020/"

####---- Load peptide data ----####

## 30min batch
f <- list.files(dataDir, pattern = "Pep.*30min.txt", 
                recursive = TRUE, full.names = TRUE)
peps30 <- read.delim(f)
## Get intensities
peps30Intensity <- peps30[, !grepl("^LFQ", colnames(peps30))]
peps30Intensity <- readSingleCellExperiment(peps30Intensity,
                                            ecol = grep("^Intensity.", colnames(peps30Intensity)),
                                            fnames = "Sequence")
colnames(peps30Intensity) <- gsub("Intensity.", "", colnames(peps30Intensity))
colnames(peps30Intensity) <- paste0("30min_", colnames(peps30Intensity))
## Get LFQ normalized intensities
peps30lfq <- peps30[, !grepl("^Intensity.\\d", colnames(peps30))]
peps30lfq <- readSingleCellExperiment(peps30lfq,
                                      ecol = grep("^LFQ", colnames(peps30lfq)),
                                      fnames = "Sequence")
colnames(peps30lfq) <- gsub("LFQ.intensity.", "", colnames(peps30lfq))
colnames(peps30lfq) <- paste0("30min_", colnames(peps30lfq))

## 60min batch
f <- list.files(dataDir, pattern = "Pep.*1hr.txt", 
                recursive = TRUE, full.names = TRUE)
peps60 <- read.delim(f)
## Get intensities
peps60Intensity <- readSingleCellExperiment(peps60,
                                            ecol = grep("^Intensity.", colnames(peps60)),
                                            fnames = "Sequence")
colnames(peps60Intensity) <- gsub("Intensity.", "", colnames(peps60Intensity))
colnames(peps60Intensity) <- paste0("60min_", colnames(peps60Intensity))

####---- Load protein data ----####

## 30min batch
f <- list.files(dataDir, pattern = "Prot.*30min.txt", 
                recursive = TRUE, full.names = TRUE)
prots30 <- read.delim(f)
## Remove unnecessary columns
prots30 <- prots30[, !grepl("^Identif|^MS.MS.count|^Peptides.\\d|nique.pep|^Sequence.coverage", colnames(prots30))]
## Get intensities
prots30Intensity <- prots30[, !grepl("^iBAQ|^LFQ", colnames(prots30))]
prots30Intensity <- readSingleCellExperiment(prots30Intensity,
                                             ecol = grep("^Intensity.", colnames(prots30Intensity)),
                                             fnames = "Protein.IDs")
colnames(prots30Intensity) <- gsub("Intensity.", "", colnames(prots30Intensity))
colnames(prots30Intensity) <- paste0("30min_", colnames(prots30Intensity))
## Get iBAQ normalized intensities
prots30ibaq <- prots30[, !grepl("^Intensity.|^LFQ", colnames(prots30))]
prots30ibaq <- readSingleCellExperiment(prots30ibaq,
                                        ecol = grep("^iBAQ.", colnames(prots30ibaq)),
                                        fnames = "Protein.IDs")
colnames(prots30ibaq) <- gsub("iBAQ.", "", colnames(prots30ibaq))
colnames(prots30ibaq) <- paste0("30min_", colnames(prots30ibaq))
## Get LFQ normalized intensities
prots30lfq <- prots30[, !grepl("^Intensity.|^iBAQ", colnames(prots30))]
prots30lfq <- readSingleCellExperiment(prots30lfq,
                                        ecol = grep("^LFQ.", colnames(prots30lfq)),
                                        fnames = "Protein.IDs")
colnames(prots30lfq) <- gsub("LFQ.intensity.", "", colnames(prots30lfq))
colnames(prots30lfq) <- paste0("30min_", colnames(prots30lfq))

## 60min batch
f <- list.files(dataDir, pattern = "Prot.*1hr.txt", 
                recursive = TRUE, full.names = TRUE)
prots60 <- read.delim(f)
## Remove unnecessary columns
prots60 <- prots60[, !grepl("^Identif|^MS.MS.count|^Peptides.\\d|nique.pep|^Sequence.coverage", colnames(prots60))]
## Get intensities
prots60Intensity <- prots60[, !grepl("^iBAQ|^LFQ", colnames(prots60))]
prots60Intensity <- readSingleCellExperiment(prots60Intensity,
                                             ecol = grep("^Intensity.", colnames(prots60Intensity)),
                                             fnames = "Protein.IDs")
colnames(prots60Intensity) <- gsub("Intensity.", "", colnames(prots60Intensity))
colnames(prots60Intensity) <- paste0("60min_", colnames(prots60Intensity))
## Get iBAQ normalized intensities
prots60ibaq <- prots60[, !grepl("^Intensity.|^LFQ", colnames(prots60))]
prots60ibaq <- readSingleCellExperiment(prots60ibaq,
                                        ecol = grep("^iBAQ.", colnames(prots60ibaq)),
                                        fnames = "Protein.IDs")
colnames(prots60ibaq) <- gsub("iBAQ.", "", colnames(prots60ibaq))
colnames(prots60ibaq) <- paste0("60min_", colnames(prots60ibaq))
## Get LFQ normalized intensities
prots60lfq <- prots60[, !grepl("^Intensity.|^iBAQ", colnames(prots60))]
prots60lfq <- readSingleCellExperiment(prots60lfq,
                                       ecol = grep("^LFQ.", colnames(prots60lfq)),
                                       fnames = "Protein.IDs")
colnames(prots60lfq) <- gsub("LFQ.intensity.", "", colnames(prots60lfq))
colnames(prots60lfq) <- paste0("60min_", colnames(prots60lfq))

####---- Create sample annotations ----####

## Sample annotation are infered from the paper and the sample names
sampleNames <- unique(c(colnames(peps30Intensity), 
                        colnames(prots30Intensity), 
                        colnames(peps60Intensity), 
                        colnames(prots60Intensity)))
annot <- DataFrame(row.names = sampleNames)
annot$LCgradient <- sub("_.*$", "", rownames(annot))
annot$NbCells <- as.numeric(sub(".*min_(\\d{1,2}).*$", "\\1", rownames(annot)))
annot$SampleType <- "MCF10A"

####---- Create the QFeatures object ----####

el <- ExperimentList(peptides_30min_intensity = peps30Intensity,
                     peptides_30min_LFQ = peps30lfq,
                     peptides_60min_intensity = peps60Intensity,
                     proteins_30min_intensity = prots30Intensity,
                     proteins_30min_iBAQ = prots30ibaq,
                     proteins_30min_LFQ = prots30lfq,
                     proteins_60min_intensity = prots60Intensity,
                     proteins_60min_iBAQ = prots60ibaq,
                     proteins_60min_LFQ = prots60lfq)
williams2020_lfq <- QFeatures(el, colData = annot)

## Create assay links
## 30 min batch
williams2020_lfq <- addAssayLink(williams2020_lfq, 
                                 from = c("peptides_30min_intensity", "peptides_30min_LFQ"), 
                                 to = "proteins_30min_intensity",
                                 varFrom = rep("Leading.razor.protein", 2), 
                                 varTo = "Majority.protein.IDs")
williams2020_lfq <- addAssayLink(williams2020_lfq, 
                                 from = c("peptides_30min_intensity", "peptides_30min_LFQ"), 
                                 to = "proteins_30min_iBAQ",
                                 varFrom = rep("Leading.razor.protein", 2), 
                                 varTo = "Majority.protein.IDs")
williams2020_lfq <- addAssayLink(williams2020_lfq, 
                                 from = c("peptides_30min_intensity", "peptides_30min_LFQ"), 
                                 to = "proteins_30min_LFQ",
                                 varFrom = rep("Leading.razor.protein", 2), 
                                 varTo = "Majority.protein.IDs")
## 60 min batch
williams2020_lfq <- addAssayLink(williams2020_lfq, 
                                 from = "peptides_60min_intensity", 
                                 to = "proteins_60min_intensity",
                                 varFrom = "Leading.razor.protein", 
                                 varTo = "Majority.protein.IDs")
williams2020_lfq <- addAssayLink(williams2020_lfq, 
                                 from = "peptides_60min_intensity", 
                                 to = "proteins_60min_iBAQ",
                                 varFrom = "Leading.razor.protein", 
                                 varTo = "Majority.protein.IDs")
williams2020_lfq <- addAssayLink(williams2020_lfq, 
                                 from = "peptides_60min_intensity", 
                                 to = "proteins_60min_LFQ",
                                 varFrom = "Leading.razor.protein", 
                                 varTo = "Majority.protein.IDs")

### Save data as Rda file
## Note: saving is assumed to occur in "scpdata/inst/scripts"
save(williams2020_lfq, 
     file = file.path("~/PhD/.localdata/scpdata/williams2020_lfq.Rda"),
     compress = "xz", 
     compression_level = 9)
