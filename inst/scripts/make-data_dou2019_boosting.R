
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
    ## Identification data (.mzid files) were downloaded at 
    ## ftp://massive.ucsd.edu/MSV000084110/result/Result_Files/
    list.files(path = dataDir, 
               pattern = paste0(batch, ".*mzid$"), 
               full.names = TRUE) %>%
        ## Read in files using `mzR`
        openIDfile %>%
        psms %>%
        ## Rename column to match with quantification data
        dplyr::rename(ScanNumber = scan.number.s.,
                      ## Avoid forbidden names for SingleCellExperiment class
                      .start = start,
                      .end = end) %>% 
        ## Remove duplicate rows
        distinct -> 
        mzid
    ## Quantification data (.txt files) were downloaded at
    ## ftp://massive.ucsd.edu/MSV000084110/other/MASIC_ReporterIons/
    list.files(path = dataDir, 
               pattern = paste0(batch, ".*txt$"), 
               full.names = TRUE) %>%
        read.table(header = TRUE, sep = "\t") %>%
        ## Add the batch name 
        mutate(Batch = batch) ->
        quant
    ## Combine data 
    left_join(mzid, quant, by = "ScanNumber")
}) %>%
    ## Bind the different batches
    bind_rows ->
    dat

## PSM is mapped to multiple isoforms, we add the isoform subscript to
## the peptide sequence

## Sample annotation provided in Table S3 in the Supplementary information
data.frame(Batch = batches,
           DatasetID = as.character(unique(dat$Dataset)),
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

## Create the `QFeatures` object
dou2019_boosting <- readSCP(dat, 
                            colDat, 
                            batchCol = "Batch", 
                            channelCol = "Channel")

####---- Peptide data ----####

## We generate the peptide data ourself, through median aggregation
dou2019_boosting <- 
    aggregateFeaturesOverAssays(dou2019_boosting, 
                                i = seq_along(dou2019_boosting),
                                fcol = "sequence", 
                                name = paste0("pep_", names(dou2019_boosting)),
                                fun = colMedians)
## We join the peptide data in a single assay
dou2019_boosting <- joinAssays(dou2019_boosting, 
                               i = grep("^pep_", names(dou2019_boosting)), 
                               name = "peptides")
## We remove the intermediate aggregation assays
dou2019_boosting <- dou2019_boosting[, , !grepl("^pep_", names(dou2019_boosting))]

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
colnames(prots) <- ifelse(grepl("Ion", channel), paste0(batch, "_", channel),
                          channel)

## Extract and format the protein expression data 
prots <- readSingleCellExperiment(prots, 
                                  ecol = grep("Ion", colname),
                                  fnames = "Protein")

## Add assay and AssayLinks to the dataset
addAssay(dou2019_boosting, prots, name = "proteins") %>%
    addAssayLink(from = "peptides", to = "proteins", 
                 varFrom = "DatabaseAccess", 
                 varTo = "Protein") ->
    dou2019_boosting

## Save data as Rda file
## Note: saving is assumed to occur in "scpdata/inst/scripts"
save(dou2019_boosting,
     compress = "xz", 
     compression_level = 9,
     file = file.path("../extdata/scpdata/dou2019_boosting.Rda"))
