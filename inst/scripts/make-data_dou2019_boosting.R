
####---- Dou et al. 2019 - testing boosting ratios ----####


## Dou, Maowei, Geremy Clair, Chia-Feng Tsai, Kerui Xu, William B. Chrisler, 
## Ryan L. Sontag, Rui Zhao, et al. 2019. “High-Throughput Single Cell Proteomics 
## Enabled by Multiplex Isobaric Labeling in a Nanodroplet Sample Preparation 
## Platform.” Analytical Chemistry, September. 
## https://doi.org/10.1021/acs.analchem.9b03349.

## This article contains 3 datasets. The script will focus on the testing of
## boosting ratios 

library(openxlsx)
library(scp)
library(mzR)
library(tidyverse)
dataDir <- "~/PhD/.localdata/SCP/dou2019_boosting/"

####---- PSM data ----####

## Load and combine the identification and quantification data
list.files(dataDir) %>%
    sub(pattern = "_[a-z]*[.][a-z]*$", 
        replacement = "", ignore.case = TRUE) %>%
    unique %>%
    grep(pattern = "xlsx", invert = TRUE, value = TRUE)->
    batches
## For every experiment (= MS run)
lapply(batches, function(batch){
    ## Identification data (.mzid files) were downloaded from
    ## ftp://massive.ucsd.edu/MSV000084110/result/Result_Files/
    ## Read in files using `mzR`
    idFile <- openIDfile(list.files(path = dataDir, 
                                    pattern = paste0(batch, ".*mzid$"), 
                                    full.names = TRUE))
    ## Extract the identification output
    matches <- psms(idFile)
    ## Extract the identification scores
    scores <- score(idFile)
    ## Combine the two tables
    matches <- left_join(matches, scores, by = "spectrumID")
    ## Rename columns
    matches <- dplyr::rename(matches,
                             ## Match the scan number column with
                             ## the quantification data
                             ScanNumber = scan.number.s.,
                             ## Avoid names forbidden by the 
                             ## SingleCellExperiment class
                             .start = start,
                             .end = end)
    ## For some unknown reason, there can be multiple E-values and 
    ## MS.GF.PepQValue assigned to a given match (and score). To avoid
    ## this, we take the lowest highest value
    matches <- group_by(matches, spectrumID, sequence) %>% 
        mutate(MS.GF.EValue = max(MS.GF.EValue),
               MS.GF.PepQValue = max(MS.GF.PepQValue))
    ## Remove duplicate rows
    matches <- distinct(matches)
    
    ## Quantification data (.txt files) were downloaded from
    ## ftp://massive.ucsd.edu/MSV000084110/other/MASIC_ReporterIons/
    list.files(path = dataDir, 
               pattern = paste0(batch, ".*txt$"), 
               full.names = TRUE) %>%
        read.table(header = TRUE, sep = "\t") %>%
        ## Add the batch name 
        mutate(Batch = batch) ->
        quant
    ## Combine data 
    left_join(matches, quant, by = "ScanNumber")
}) %>%
    ## Bind the different batches
    bind_rows ->
    psms
## Include also contaminant information
psms$isContaminant <- grepl("Contaminant", psms$DatabaseAccess)

## PSM data is mapped to multiple peptides and peptides to multiple 
## proteins. In such cases, the PSM quantitative data is duplicated.
## To avoid this duplication, we combine all features belonging to the
## same spectra. This means we are creating peptide groups when a 
## spectrum is matched to multiple peptides, and proteins groups when 
## a peptide is mapped to multiple proteins. 

## When combining multiple lines into a single PSM, some values differ.
## If so, we concatenate those values in a single entry. Note that when
## one of the peptides is a decoy, we flag the corresponding spectra 
## as a decoy. 
combineFeatures <- function(x) {
    if (length(x) == 1) return(x)
    ## isDecoy is a logical. if any PSM is a decoy, the spectrum is 
    ## considered as decoy
    if (is.logical(x)) return(all(x)) 
    xuni <- unique(x)
    if (length(xuni) == 1) return(xuni)
    ## When a spectrum matches to multiples proteins, we combine those
    ## in a protein group
    paste0(xuni, collapse = ";")
}
## In order to concatenate some columns that are not unique for a 
## given spectra, we need to convert those columns to characters
psms$.start <- as.character(psms$.start)
psms$.end <- as.character(psms$.end)
psms$DBseqLength <- as.character(psms$DBseqLength)
psms$modNum <- as.character(psms$modNum)
psms$calculatedMassToCharge <- as.character(psms$calculatedMassToCharge)
## Combine the features as single spectrum information
psms <- group_by(psms, spectrumID, Batch)
psms <- dplyr::summarise(psms, across(.fn = combineFeatures))

####---- Sample annotation ----####

## Sample annotation provided in Table S3 in the Supplementary information
data.frame(Batch = batches,
           DatasetID = as.character(unique(psms$Dataset)),
           Ion_126.128 = c("Raw", "SVEC", "C10", "SVEC", "C10", "SVEC"),
           Ion_127.125 = c("Raw", "SVEC", "C10", "SVEC", "C10", "SVEC"),
           Ion_127.131 = c("Raw", "SVEC", "Raw", "SVEC", "C10", "SVEC"),
           Ion_128.128 = c("SVEC", "C10", "Raw", "C10", "Raw", "C10"),
           Ion_128.134 = c("SVEC", "C10", "Raw", "C10", "Raw", "C10"),
           Ion_129.131 = c("SVEC", "C10", "SVEC", "C10", "Raw", "C10"),
           Ion_129.138 = c("C10", "Raw", "SVEC", "Raw", "SVEC", "Raw"),
           Ion_130.135 = c("C10", "Raw", "empty", "empty", "empty", "empty"),
           Ion_130.141 = c("C10", "Raw", "C10", "Raw", "SVEC", "Raw"),
           Ion_131.138 = c("empty", "empty", "boost", "boost", "boost", "boost")) %>%
    pivot_longer(cols = matches("Ion"), 
                 names_to = "Channel", 
                 values_to = "SampleType") ->
    colDat

####---- Create the QFeatures object ----####

## Create the `QFeatures` object
dou2019_boosting <- readSCP(psms, 
                            colDat, 
                            batchCol = "Batch", 
                            channelCol = "Channel")

####---- Peptide data ----####

## We generate the peptide data ourself

## First, we filter out low quality PSM
dou2019_filt <- filterFeatures(dou2019_boosting, ~ 
                                   ! isContaminant & 
                                   ! isDecoy & 
                                   MS.GF.QValue < 0.01)
## Peptide data is computed by median aggregation of PSMs to peptides
dou2019_filt <- 
    aggregateFeaturesOverAssays(dou2019_filt, 
                                i = seq_along(dou2019_filt),
                                fcol = "sequence", 
                                name = paste0("pep_", names(dou2019_filt)),
                                fun = colMedians)
## We join the peptide data in a single assay
dou2019_filt <- joinAssays(dou2019_filt, 
                           i = grep("^pep_", names(dou2019_filt)), 
                           name = "peptides")
## We remove the peptide groups
sce <- dou2019_filt[["peptides"]]
sce <- sce[!grepl("[;]", rownames(sce)), ]
## We add the assay to the exported dataset and create the links 
## between PSM and peptide data 
dou2019_boosting <- addAssay(dou2019_boosting, sce, name = "peptides")
dou2019_boosting <- 
    addAssayLink(dou2019_boosting, from = 1:6, to = "peptides", 
                 varFrom = rep("sequence", 6), varTo = "sequence")

####---- Protein data ----####

## Load the data
## Data was downloaded from https://doi.org/10.1021/acs.analchem.9b03349.
datXlsx <- loadWorkbook(list.files(path = dataDir, 
                                   pattern = "xlsx$",  
                                   full.names = TRUE))
sheets <- c(2, ## 0ng boosting
            4, ## 5ng boosting
            6) ## 50ng boosting
lapply(sheets, 
       read.xlsx,
       xlsxFile = datXlsx, colNames = TRUE, startRow = 2) %>%
    ## Remove the peptide and spectra counts
    lapply(function(x) x[, -(2:3)]) %>%
    ## Joint the datasets based on the protein
    reduce(full_join, by = "Protein") ->
    prots

## Rename columns
channel <- sub(pattern = "_Data.*$", replacement = "", colnames(prots))
datasetID <- sub(pattern = "Ion.*ID_", replacement = "",  colnames(prots))
batch <- colDat$Batch[match(datasetID, colDat$DatasetID)]
colnames(prots) <- ifelse(grepl("Ion", channel), paste0(batch, channel),
                          channel)

## Extract and format the protein expression data 
prots <- readSingleCellExperiment(prots, 
                                  ecol = grep("Ion", colnames(prots)),
                                  fnames = "Protein")

## Add assay and AssayLinks to the dataset
dou2019_boosting <- addAssay(dou2019_boosting, prots, name = "proteins")
dou2019_boosting <- addAssayLink(dou2019_boosting,
                                 from = "peptides", to = "proteins", 
                                 varFrom = "DatabaseAccess", 
                                 varTo = "Protein")

## Save data as Rda file
## Note: saving is assumed to occur in "scpdata/inst/scripts"
save(dou2019_boosting,
     compress = "xz", 
     compression_level = 9,
     file = file.path("~/PhD/.localdata/scpdata/dou2019_boosting.Rda"))
