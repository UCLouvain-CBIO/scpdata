
####---- Zhu et al. 2018, Nature Communications - HeLa dilutions ----####

## Zhu, Ying, Paul D. Piehowski, Rui Zhao, Jing Chen, Yufeng Shen, Ronald J. 
## Moore, Anil K. Shukla, et al. 2018. “Nanodroplet Processing Platform for Deep 
## and Quantitative Proteome Profiling of 10-100 Mammalian Cells.” Nature 
## Communications 9 (1): 882.

library(SingleCellExperiment)
library(scp)
library(tidyverse)
dataDir <- "~/PhD/.localdata/SCP/zhu2018NC_hela/"

## The data was downloaded from ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2018/01/PXD006847
## to scpdata/inst/extdata/zhu2018NC_hela

####---- Get the annotation data ----####

## The annotation file was manually created from the samples names and
## the available information from the methods section
annot <- read.csv(paste0(dataDir, "sample_annotation.csv"), 
                  row.names = "SampleName")

####---- Peptide data ----####

## Load the quantification data
peps <- read.table(paste0(dataDir, "CulturedCells_peptides.txt"),
                   sep = "\t", header = TRUE)
peps <- readSingleCellExperiment(peps,
                                 ecol = grep("^Intensity.", colnames(peps)),
                                 fnames = "Sequence")
colnames(peps) <- gsub("Intensity.", "", colnames(peps))

####---- Protein data ----####

## Load the quantification data
prots <- read.table(paste0(dataDir, "CulturedCells_proteinGroups.txt"),
                    sep = "\t", header = TRUE)
## Remove unnecessary columns
sel <- !grepl("Peptides.*[HMT]|^Identif|^Razor.*[HMT]|Sequence.*[HMT]|Unique.*[HMT]|MS.MS.*[HMT]",
              colnames(prots))
prots <- prots[, sel]
## Split protein data based on the quantification method:
## 1. Protein intensity
protsInt <- prots[, !grepl("^iBAQ|^LFQ", colnames(prots))]
protsInt <- readSingleCellExperiment(protsInt, 
                                     ecol = grep("^Intensity.", colnames(protsInt)),
                                     fnames = "Protein.IDs")
colnames(protsInt) <- gsub("Intensity.", "", colnames(protsInt))
## 2. LFQ
protsLFQ <- prots[, !grepl("^iBAQ|^Intensity", colnames(prots))]
protsLFQ <- readSingleCellExperiment(protsLFQ, 
                                     ecol = grep("^LFQ.", colnames(protsLFQ)),
                                     fnames = "Protein.IDs")
colnames(protsLFQ) <- gsub("LFQ.intensity.", "", colnames(protsLFQ))
## 3. iBAQ
protsIBAQ <- prots[, !grepl("^LFQ|^Intensity", colnames(prots))]
protsIBAQ <- readSingleCellExperiment(protsIBAQ, 
                                      ecol = grep("^iBAQ.", colnames(protsIBAQ)),
                                      fnames = "Protein.IDs")
colnames(protsIBAQ) <- gsub("iBAQ.", "", colnames(protsIBAQ))


####---- Create the QFeatures object ----####

el <- ExperimentList(peptides = peps,
                     proteins_intensity = protsInt,
                     proteins_LFQ = protsLFQ,
                     proteins_iBAQ = protsIBAQ)
zhu2018NC_hela <- QFeatures(el, colData = annot)

## Create assay links
zhu2018NC_hela <- addAssayLink(zhu2018NC_hela, 
                                 from = "peptides", to = "proteins_intensity",
                                 varFrom = "Leading.razor.protein", 
                                 varTo = "Majority.protein.IDs")
zhu2018NC_hela <- addAssayLink(zhu2018NC_hela, 
                                 from = "peptides", to = "proteins_LFQ",
                                 varFrom = "Leading.razor.protein", 
                                 varTo = "Majority.protein.IDs")
zhu2018NC_hela <- addAssayLink(zhu2018NC_hela, 
                                 from = "peptides", to = "proteins_iBAQ",
                                 varFrom = "Leading.razor.protein", 
                                 varTo = "Majority.protein.IDs")

# Save data as Rda file
# Note: saving is assumed to occur in "(...)/scpdata/inst/scripts"
save(zhu2018NC_hela, 
     file = "~/PhD/.localdata/scpdata/zhu2018NC_hela.Rda",
     compress = "xz", 
     compression_level = 9)
