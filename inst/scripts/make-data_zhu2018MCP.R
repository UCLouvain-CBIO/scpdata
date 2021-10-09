
####---- Zhu et al. 2018, Molecular & Cellular Proteomics ----####

# Zhu, Ying, Maowei Dou, Paul D. Piehowski, Yiran Liang, Fangjun Wang, Rosalie 
# K. Chu, William B. Chrisler, et al. 2018. “Spatially Resolved Proteome Mapping 
# of Laser Capture Microdissected Tissue with Automated Sample Transfer to 
# Nanodroplets.” Molecular & Cellular Proteomics: MCP 17 (9): 1864–74.

library(scp)
library(mzR)
library(tidyverse)
dataDir <- "~/PhD/.localdata/SCP/zhu2018MCP/"

# The data was downloaded from ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2018/07/PXD008844

####---- Peptide data ----####

## Load the peptide data
peps <- read.table(paste0(dataDir, "MaxQuant_Peptides.txt"),
                   sep = "\t", header = TRUE)
peps <- readSingleCellExperiment(peps,
                                 ecol = grep("^Intensity.", colnames(peps)),
                                 fnames = "Sequence")
colnames(peps) <- gsub("Intensity.", "", colnames(peps))

####---- Sample annotations ----####

## Create the sample annotation table
annot <- DataFrame(row.names = colnames(peps))
annot$SectionWidth <- sub("^([0-9]*um)_.*$", "\\1", rownames(annot))
annot$SampleType <- str_extract(rownames(annot), "CC|CP|Mix")
annot$SampleType <- ifelse(is.na(annot$SampleType), "CTX", annot$SampleType)
annot$Prep <- str_match(rownames(annot), "08.*Prep$")[, 1]
annot$Replicate <- str_match(rownames(annot), "(_|Mix)(\\d)(_|$)")[, 3]

####---- Protein data ----####

## Load the quantification data
prots <- read.table(paste0(dataDir, "MaxQuant_ProteinGroups.txt"),
                    sep = "\t", header = TRUE)
## Remove unnecessary columns
sel <- !grepl("Peptides.*\\d|^Identif|^Razor.*\\d|Sequence.*\\d|Unique.*\\d|MS.MS.*\\d",
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
zhu2018MCP <- QFeatures(el, colData = annot)

## Create assay links
zhu2018MCP <- addAssayLink(zhu2018MCP, 
                                  from = "peptides", to = "proteins_intensity",
                                  varFrom = "Leading.razor.protein", 
                                  varTo = "Majority.protein.IDs")
zhu2018MCP <- addAssayLink(zhu2018MCP, 
                                  from = "peptides", to = "proteins_LFQ",
                                  varFrom = "Leading.razor.protein", 
                                  varTo = "Majority.protein.IDs")
zhu2018MCP <- addAssayLink(zhu2018MCP, 
                                  from = "peptides", to = "proteins_iBAQ",
                                  varFrom = "Leading.razor.protein", 
                                  varTo = "Majority.protein.IDs")

### Save data as Rda file
## Note: saving is assumed to occur in "scpdata/inst/scripts"
save(zhu2018MCP, 
     file = file.path("~/PhD/.localdata/scpdata/zhu2018MCP.Rda"),
     compress = "xz", 
     compression_level = 9)
