
####---- Krull et al, 2024 ---####

## Krull, K. K., Ali, S. A., & Krijgsveld, J. 2024. "Enhanced feature matching 
## in single-cell proteomics characterizes IFN-γ response and co-existence of 
## Cell States." Nature Communications, 15(1). https://doi.org/10.1038/s41467-024-52605-x 

library(SingleCellExperiment)
library(scp)
library(tidyverse)

root <- "~/localdata/SCP/krull2024/"

## Data downloaded from:
## https://www.ebi.ac.uk/pride/archive/projects/PXD053464, 
## look for 03_SingleCell_Searches.zip. The relevant file is 
## 03_SingleCell_Searches/10cell_DIA-ME/report.tsv

report <- read.delim(paste0(root, "report.tsv"))

## Correct the colname 
report <- report %>% 
  rename(MS1.Area = Ms1.Area)

## Build the annotation table from the run names
annot <- DataFrame(File.Name = unique(report$File.Name))
annot$Run <- gsub("^D.*[\\]|[.]d$", "", annot$File.Name)
otherVars <- strsplit(annot$Run, "_")
otherVars <- lapply(otherVars, function(x) {
  if (length(x) == 12) x <- x[-4]
  x
})
otherVars <- do.call(rbind, otherVars)
colnames(otherVars) <- c("Date", "..undetermined..", "SampleAnnotation", 
                         "SampleType", "..undetermined..", "..undetermined..", 
                         "..undetermined..", "RunID")
annot <- cbind(annot, otherVars)
annot <- annot[, !grepl("\\.\\.undetermined\\.\\.", colnames(annot))]

## Format to a QFeatures object
krull2024 <- readSCPfromDIANN(annot, report)

####---- Add the protein data ----####
## Data downloaded from:
## Supplementary Data 8 of the paper 
## https://www.nature.com/articles/s41467-024-52605-x#Sec21
## 41467_2024_52605_MOESM12_ESM.xlsx and DIA-ME sheet exported as csv.

prot_data <- read.csv(paste0(root, "41467_2024_52605_MOESM12_ESM_DIA-ME.csv"))

## Create the SingleCellExperiment object
prots <- readSingleCellExperiment(prot_data, 
                                ecol = 3:145)

## Name rows with peptide sequence
rownames(prots) <- prot_data$Protein.Group

from <- names(krull2024)
krull2024 <- addAssay(krull2024, prots, "proteins")
krull2024 <- addAssayLink(krull2024, 
                            from = from, 
                            to = "proteins", 
                            varFrom = rep("Protein.Group", length(from)),
                            varTo = "Protein.Group")

## Save data
save(krull2024,
     file = file.path(paste0(root, "krull2024.Rda")),
     compress = "xz",
     compression_level = 9)

