
####---- Dou et al. 2019 - HeLa digests ----####


## Dou, Maowei, Geremy Clair, Chia-Feng Tsai, Kerui Xu, William B. Chrisler, 
## Ryan L. Sontag, Rui Zhao, et al. 2019. “High-Throughput Single Cell Proteomics 
## Enabled by Multiplex Isobaric Labeling in a Nanodroplet Sample Preparation 
## Platform.” Analytical Chemistry, September. 
## https://doi.org/10.1021/acs.analchem.9b03349.

## This article contains 3 datasets. The script will focus on the HeLa digests

library(openxlsx)
library(scp)
library(mzR)
library(tidyverse)
setwd("inst/scripts/")


####---- PSM data ----####


## Load and combine the identification and quantification data
batches <- c("Hela_run_1", "Hela_run_2") 
## For every experiment (= MS run)
lapply(batches, function(batch){
  ## Identification data (.mzid files) were downloaded at 
  ## ftp://massive.ucsd.edu/MSV000084110/result/Result_Files/
  list.files("../extdata/dou2019_lysates/", 
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
  list.files("../extdata/dou2019_lysates", 
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

## Create the sample metadata
## The table is manually created because the information is taken from the 
## supplementary information file (Figure S2)
channels <- grep("^Ion_.*\\d$", colnames(dat), value = TRUE)
colDat <- data.frame(Channel = rep(channels, length(batches)),
                     Batch = rep(batches, each = length(channels)),
                     SampleType = rep(c(rep("Lysate", 7), 
                                        rep("Blank", 2), 
                                        "Carrier"),
                                      length(batches)),
                     InputAmount_ng = rep(c(rep(0.2, 7), 
                                            rep(0, 2), 
                                            10),
                                          length(batches)))

## Create the `QFeatures` object
dou2019_lysates <- readSCP(dat, 
                        colDat, 
                        batchCol = "Batch", 
                        channelCol = "Channel")


####---- Protein data ----####


## Load the data
## Data was downloaded from https://doi.org/10.1021/acs.analchem.9b03349.
"../extdata/dou2019_lysates/ac9b03349_si_003.xlsx" %>%
  loadWorkbook %>%
  read.xlsx(sheet = 6, colNames = FALSE) -> 
  dat

## Get the column names 
## The columns names should match the column names contained in the 
## 'dou2019_lysates' object
dat[1:3, ] %>%
  t %>% 
  as_tibble %>%
  mutate(colname = sub(pattern = "_Dat.*$", replacement = "", `3`),
         batch = sub("run", "", `1`), 
         colname = ifelse(grepl("Ion", colname),
                          paste0("Hela_run_", batch, "_", colname),
                          colname)) %>%
  pull(colname) ->
  colname
## Extract and format the protein expression data 
dat[-(1:3), ] %>%
  magrittr::set_colnames(colname) %>%
  mutate_at(vars(contains("Ion")), as.numeric) %>% 
  readSingleCellExperiment(ecol = grep("Ion", colname),
                           fnames = "Protein") ->
  dat

## Add assay and AssayLinks to the dataset
addAssay(dou2019_lysates, dat, name = "proteins") %>%
  addAssayLink(from = names(dou2019_lysates)[1:2], to = "proteins", 
               varFrom = rep("DatabaseAccess", 2), 
               varTo = "Protein") ->
  dou2019_lysates

## Save data as Rda file
## Note: saving is assumed to occur in "scpdata/inst/scripts"
save(dou2019_lysates, 
     compress = "xz", 
     compression_level = 9,
     file = file.path("../EHdata/scpdata/dou2019_lysates.Rda"))

