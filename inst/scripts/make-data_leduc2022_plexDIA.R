####---- Leduc et al. 2022 ---####


## Leduc, Andrew, R. Gray Huffman, and Nikolai Slavov. 2021. “Droplet 
## Sample Preparation for Single-Cell Proteomics Applied to the Cell 
## Cycle.” bioRxiv. https://doi.org/10.1101/2021.04.24.441211.

## This is the script for generating the plexDIA dataset

library("scp")
library("SCP.replication")
library("tidyverse")

## All files were downloaded from 
## https://drive.google.com/drive/folders/117ZUG5aFIJt0vrqIxpKXQJorNtekO-BV

datadir <- "~/PhD/.localdata/SCP/leduc2022/plexDIA/"

####---- Prepare sample annotations ----####

# The sample annotations are provided in 2 separate tables:
design <- read.csv(paste0(datadir, "annotation_plexDIA.csv"))
design <- rename(design, 
                 Run = Set)

## Failed injections not to include when doing DIA-NN max LFQ
runs <-c('wAL_plexMel65', 'wAL_plexMel66', 'wAL_plexMel02',
         'wAL_plexMel04', 'wAL_plexMel06', 'wAL_plexMel08',
         'wAL_plexMel10', 'wAL_plexMel12', 'wAL_plexMel14',
         'wAL_plexMel16', 'wAL_plexMel18', 'wAL_plexMel20',
         'wAL_plexMel22', 'wAL_plexMel24', 'wAL_plexMel26',
         'wAL_plexMel28', 'wAL_plexMel30', 'wAL_plexMel32')
design <- design[!design$Run %in% runs, ]

# Clean the sample metadata so that it meets the requirements for
# `scp::readSCP`. We first need to transform the design (run x 
# reporter ion) to a long table so that one line is one sample. 
design <- pivot_longer(design, -Run, names_to = "Label", 
                       values_to = "SampleAnnotation")
design$SampleType <- recode(design$SampleAnnotation, 
                            neg = "NegativeControl",
                            M = "Melanoma")
design$Label <- sub("^X", "", design$Label)
design$File.Name <- make.names(paste0("D:\\AL\\AL\\", design$Run, ".raw"))

####---- Prepare precursor data ----####

extracted <- read.delim(paste0(datadir, "report.pr_matrix_channels_ms1_extracted.tsv"))
extracted$Stripped.Sequence.Charge <- paste0(extracted$Stripped.Sequence, extracted$Precursor.Charge)
report <- read.delim(paste0(datadir, "report_plexDIA_mel_nPOP.tsv"))
report$File.Name <- make.names(report$File.Name)
leduc2022_plexDIA <- readSCPfromDIANN(colData = design, 
                                      reportData = report,
                                      extractedData = extracted,
                                      multiplexing = "mTRAQ")

####---- Add peptide data ----####

peps <- read.csv(paste0(datadir, "plexDIA_peptide.csv"))
peps <- rename(peps, Stripped.Sequence.Charge = X)
peps <- readSingleCellExperiment(peps, fname = "Stripped.Sequence.Charge",
                                 ecol = grep("^wAL", colnames(peps)))
colnames(peps) <- make.names(paste0("D:\\AL\\AL\\", colnames(peps)))
leduc2022_plexDIA <- addAssay(leduc2022_plexDIA, peps, name = "peptides")
leduc2022_plexDIA <- addAssayLink(leduc2022_plexDIA, 
                                  from = "Ms1Extracted", to = "peptides",
                                  varFrom = "Stripped.Sequence.Charge",
                                  varTo = "Stripped.Sequence.Charge")


####---- Add protein data ----####

prots <- read.csv(paste0(datadir, "plexDIA_protein_imputed.csv"))
prots$Protein.Group <- rownames(prots)
prots <- readSingleCellExperiment(prots, fname = "Protein.Group",
                                 ecol = grep("^wAL", colnames(prots)))
colnames(prots) <- make.names(paste0("D:\\AL\\AL\\", colnames(prots)))
leduc2022_plexDIA <- addAssay(leduc2022_plexDIA, prots, name = "proteins")
leduc2022_plexDIA <- addAssayLink(leduc2022_plexDIA, 
                                  from = "Ms1Extracted", to = "proteins",
                                  varFrom = "Protein.Group",
                                  varTo = "Protein.Group")

## Save data as Rda file
save(leduc2022_plexDIA, 
     file = "~/PhD/.localdata/scpdata/leduc2022_plexDIA.Rda",
     compress = "xz", 
     compression_level = 9)
