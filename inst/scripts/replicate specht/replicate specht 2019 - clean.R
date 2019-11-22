
####---- Description ----#### 

# This script will replicate the analysis performed in Specht et al. 2019. The 
# processing of the MaxQuant output file will be reproduced but using clean  
# function implemented in the framework of MSnbase. 


####---- Load the data ----####

library(MSnbase)
setwd("./inst/scripts/")

# Load meta data
samp <- read.csv("../extdata/specht2019/annotation_fp60-97.csv", row.names = 1)
batch <- read.csv("../extdata/specht2019/batch_fp60-97.csv", row.names = 1)

# Load the data and create an MSnSet
mq_file <- "../extdata/specht2019/ev_updated.txt"
coln <- colnames(read.table(mq_file, header = TRUE, sep = "\t", nrow = 1))
sc <- readMSnSet2(file = mq_file, fnames = "id", sep = "\t", header = TRUE,
                  ecol = grep("intensity[.]\\d", coln, value = TRUE))

# Keep only sc runs 
sel <- !grepl("blank|_QC_|col[12][0-9]", fData(sc)$Raw.file)
sc <- sc[sel, ]


####---- Filter identification ----####

sc %>% filter(TRUE)

# Make sure all runs are described in design, if not, print and remove them:
not.described<-unique(ev$Raw.file)[ !unique(ev$Raw.file) %in% colnames(design) ]
ev<-ev[!ev$Raw.file%in%not.described,]

# Filter out reverse hits, contaminants, and contaminated spectra... 
ev<-ev[which(ev$Reverse!="+"),]
ev<-ev[-grep("REV", ev$Leading.razor.protein),]
ev<-ev[ev$Potential.contaminant!="+",]
ev<-ev[ev$PIF>0.8,]
ev<-ev[!is.na(ev$PIF),]



