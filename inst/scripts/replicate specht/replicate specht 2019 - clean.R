
####---- Description ----#### 

# This script will replicate the analysis performed in Specht et al. 2019. The 
# processing of the MaxQuant output file will be reproduced but using clean  
# function implemented in the framework of MSnbase. 


####---- Load the data ----####

library(MSnbase)
setwd("./inst/scripts/")

# Load meta data
samp <- read.csv("../extdata/specht2019/annotation_fp60-97.csv", 
                 row.names = 1, check.names = F)
batch <- read.csv("../extdata/specht2019/batch_fp60-97.csv", row.names = 1)

# Load the data and create an MSnSet
mq_file <- "../extdata/specht2019/ev_updated.txt"
coln <- colnames(read.table(mq_file, header = TRUE, sep = "\t", nrow = 1))
sc <- sc0 <- readMSnSet2(file = mq_file, fnames = "id", sep = "\t", header = TRUE,
                         ecol = grep("intensity[.]\\d", coln, value = TRUE))

# Keep only sc runs and runs that have associated sample metadata
sel <- !grepl("blank|_QC_|col[12][0-9]", fData(sc)$Raw.file)
sel <- sel & fData(sc)$Raw.file %in% colnames(samp)
sc <- sc[sel, ]


####---- Filter based on identification measures ----####


# Remove the reverse hits (from decoy database) and contaminants
sc <- sc[fData(sc)$Reverse != "+", ]
sc <- sc[!grepl("^REV", fData(sc)$Leading.razor.protein), ]
sc <- sc[fData(sc)$Potential.contaminant != "+", ]
sc <- sc[fData(sc)$PIF > 0.8 & !is.na(fData(sc)$PIF), ]

# Remove spectra with poor identification confidence
# The PEP and q-values were updated using DART-ID
# TODO discuss with Laurent ! Why not deleting peptides with high FDR instead of deleting the proteins for which all associated peptides have high FDR 
qprots <- unique(fData(sc)$Leading.razor.protein[fData(sc)$dart_qval < 0.01])
sc <- sc[fData(sc)$Leading.razor.protein %in% qprots, ]
sc <- sc[fData(sc)$dart_PEP < 0.02, ]

# Remove runs with insufficient identified peptides 
pep.t <- table(fData(sc)$Raw.file)
sc <- sc[fData(sc)$Raw.file %in% names(pep.t[pep.t >= 300]),] # 25 runs were removed


####--- Filter based on sample to carrier ratio


