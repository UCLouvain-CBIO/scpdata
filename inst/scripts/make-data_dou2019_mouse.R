
####---- Dou et al. 2019 - mouse cell data set ----####


## Dou, Maowei, Geremy Clair, Chia-Feng Tsai, Kerui Xu, William B. Chrisler, 
## Ryan L. Sontag, Rui Zhao, et al. 2019. “High-Throughput Single Cell Proteomics 
## Enabled by Multiplex Isobaric Labeling in a Nanodroplet Sample Preparation 
## Platform.” Analytical Chemistry, September. 
## https://doi.org/10.1021/acs.analchem.9b03349.

## This article contains 3 datasets. The script will focus on the profiling of 
## murine cell populations

library(openxlsx)
library(scp)
library(mzR)
library(tidyverse)
setwd("inst/scripts/")
dataDir <- "../extdata/dou2019_mouse/"

####---- PSM data ----####


## Load and combine the identification and quantification data
list.files(dataDir) %>%
  sub(pattern = "_[a-z]*[.][a-z]*$", 
      replacement = "", ignore.case = TRUE) %>%
  unique ->
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

## Sample annotation provided in Table S3 in the Supplementary information
data.frame(Batch = batches,
           DatasetID = as.character(unique(dat$Dataset)),
           Ion_126.128 = c(rep("C10", 6), rep("SVEC", 6)),
           Ion_127.125 = c(rep("C10", 3), rep("SVEC", 6), rep("RAW", 3)),
           Ion_127.131 = c(rep("SVEC", 6), rep("RAW", 6)),
           Ion_128.128 = c(rep("SVEC", 3), rep("RAW", 6), rep("C10", 3)),
           Ion_128.134 = c(rep("RAW", 6), rep("C10", 6)),
           Ion_129.131 = c(rep("RAW", 3), rep("C10", 6), rep("SVEC", 3)),
           Ion_129.138 = c(rep("reference", 12)),
           Ion_130.135 = c(rep("empty", 12)),
           Ion_130.141 = c(rep("empty", 12)),
           Ion_131.138 = c(rep("boost", 12)), 
           ## Note there is a mistake in the table. The author provide 
           ## the annotation for 10 channles whereas 11 channels are recorded.
           ## The last column in the data has little signal, so sample type is 
           ## assumed to be "empty"
           Ion_131.144 = c(rep("empty", 12))) %>%
  pivot_longer(cols = matches("Ion"), 
               names_to = "Channel", 
               values_to = "SampleType") ->
  colDat

## Create the `QFeatures` object
dou2019_mouse <- readSCP(dat, 
                         colDat, 
                         batchCol = "Batch", 
                         channelCol = "Channel")


####---- Protein data ----####


## Load the data
## Data was downloaded from https://doi.org/10.1021/acs.analchem.9b03349.
list.files(path = dataDir, 
           pattern = "xlsx$",  
           full.names = TRUE) %>%
  loadWorkbook %>%
  read.xlsx(sheet = 2, colNames = FALSE) ->
  dat

## Get the column names 
## The columns names should match the column names contained in the 
## 'dou2019_hela' object
dat[1:2, ] %>%
  t %>% 
  as_tibble %>%
  mutate(DatasetID = sub(pattern = "Ion_.*ID_", replacement = "", `2`),
         Batch = colDat$Batch[match(DatasetID, colDat$DatasetID)],
         colname = sub(pattern = "_Dat.*$", replacement = "", `2`),
         colname = ifelse(grepl("Ion", colname),
                          paste0(Batch, "_", colname),
                          colname)) %>%
  pull(colname) ->
  colname

## Extract and format the protein expression data 
dat[-(1:2), ] %>%
  set_colnames(colname) %>%
  mutate_at(vars(contains("Ion")), as.numeric) %>% 
  readSingleCellExperiment(ecol = grep("Ion", colname),
                           fnames = "Protein") ->
  dat

## Add assay and AssayLinks to the dataset
addAssay(dou2019_mouse, dat, name = "proteins") %>%
  addAssayLink(from = 1:12, to = "proteins", 
               varFrom = rep("DatabaseAccess", 12), 
               varTo = "Protein") ->
  dou2019_mouse

## Save data as Rda file
## Note: saving is assumed to occur in "scpdata/inst/scripts"
save(dou2019_mouse,
     compress = "xz", 
     compression_level = 9,
     file = file.path("../EHdata/dou2019_mouse.rda"))
