
####---- Schoof et al. 2021 ---####


## Schoof, Erwin M., Benjamin Furtwängler, Nil Üresin, Nicolas Rapin,
## Simonas Savickas, Coline Gentil, Eric Lechman, Ulrich auf Dem 
## Keller, John E. Dick, and Bo T. Porse. 2021. “Quantitative 
## Single-Cell Proteomics as a Tool to Characterize Cellular 
## Hierarchies.” Nature Communications 12 (1): 745679.

library(SingleCellExperiment)
library(scp)
library(tidyverse)

## The data files were retrieved from the `SCeptre_FINAL.zip` file 
## provided in the PRIDE repository PXD020586
## https://www.ebi.ac.uk/pride/archive/projects/PXD020586
## We unziped the file and retrieved the data folder located at 
## `SCeptre/Schoof_et_al/data` (here renamed`.../schoof2021/`)

dataDir <- "~/PhD/.localdata/SCP/schoof2021/bulk/"
dataDir <- "~/PhD/tmp/SCeptre/Schoof_et_al/data/bulk/"
## Load protein data (to get sample names)
list.files(path = dataDir, pattern = "Proteins.txt", full.names = TRUE) %>% 
    read.table(header = TRUE, sep = "\t") ->
    prot


####---- Create sample annotation ----####

## Extract annotations from the sample names
data.frame(Colname = colnames(prot)) %>%
    filter(grepl("^Abundances",  Colname)) %>%
    mutate(File.ID = sub("^.*(F[0-9]*)[.].*$", "\\1", Colname),
           Channel = sub("^.*[.](.*)$", "\\1", Colname)) %>%
    select(-Colname) ->
    meta

## Get the ProteomeDiscoverer metadata
list.files(path = dataDir, pattern = "InputFiles.txt", full.names = TRUE) %>%
    read.table(header = TRUE, 
               sep = "\t") %>%
    mutate(File.Name = sub("^.*[\\](.*)$",  "\\1",  File.Name)) ->
    pdData
## Add to the metadata table
meta <- left_join(meta, pdData, by = "File.ID")

## Get the file to plate mapping
list.files(path = dataDir, 
           pattern = "file_sample",
           full.names = TRUE) %>% 
    read.table(header = TRUE, sep = "\t") ->
    fileMap
## Add to the metadata table
meta <- left_join(meta, fileMap, by = "File.Name")

## Get the plate to layouts mappings
list.files(path = dataDir, pattern = "plate_layout", full.names = TRUE) %>% 
    read.table(header = TRUE, sep = "\t") ->
    plateLayouts

## Get the sample labelling layouts
lapply(list.files(path = dataDir, pattern = "sample_layout"), 
       function(file) {
           df <- read.table(paste0(dataDir, file), header = TRUE, sep = "\t",
                            check.names = FALSE)
           colnames(df)[1] <- "Row"
           df <- pivot_longer(df, -Row, names_to = "Col", values_to = "Sample")
           df$Sample.Layout <- file
           df
       }) %>% 
    do.call(what = rbind) ->
    sampleLayouts
## Add to the plate data
plateLayouts <- merge(plateLayouts, sampleLayouts, by = "Sample.Layout")

## Get the sample sorting layouts
lapply(list.files(path = dataDir, pattern = "sort_layout"), 
       function(file) {
           df <- read.table(paste0(dataDir, file), header = TRUE, sep = "\t",
                            check.names = FALSE)
           colnames(df)[1] <- "Row"
           df <- pivot_longer(df, -Row, names_to = "Col", values_to = "SortedPopulation")
           df$Sort.Layout <- file
           df
       }) %>% 
    do.call(what = rbind) ->
    sortLayouts
## Add to the plate data
plateLayouts <- merge(plateLayouts, sortLayouts, by = c("Row", "Col", "Sort.Layout"))

## Get the sample labelling layouts
lapply(list.files(path = dataDir, pattern = "label_layout"), 
       function(file) {
           df <- read.table(paste0(dataDir, file), header = TRUE, sep = "\t",
                            check.names = FALSE)
           colnames(df)[1] <- "Row"
           df <- pivot_longer(df, -Row, names_to = "Col", values_to = "Channel")
           df$Channel[df$Channel == ""] <- NA
           df$Label.Layout <- file
           df
       }) %>% 
    do.call(what = rbind) ->
    labelLayouts
## Add to the plate data
plateLayouts <- merge(plateLayouts, labelLayouts, by = c("Row", "Col", "Label.Layout"))

## Get the FACS data and combine in a single table
lapply(list.files(path = dataDir, pattern = "facs"), 
       function(file) {
           df <- read.table(paste0(dataDir, file), header = TRUE, sep = "\t")   
           df$Facs.Data <- file
           df
       }) %>% 
    do.call(what = rbind) %>%
    select(-Row, -Column) %>% 
    dplyr::rename(CD34_APC.Cy7.A = APC.Cy7.A,
                  CD38_PE.A = PE.A) ->
    facsData
## Add to the plate data
plateLayouts$Well <- paste0(plateLayouts$Row, plateLayouts$Col)
plateLayouts <- merge(plateLayouts, facsData, by = c("Well", "Facs.Data"))

## Add all the layouts to the metadata table
## Note we here use left_join because some Sample are NA, but we still
## want to keep them. 
meta <- left_join(meta, plateLayouts, by = c("Channel", "Plate", "Sample"))

## Add further info retrieved from the article
meta$SampleType <- "empty"
meta$SampleType[meta$Channel == "126" & grepl("bulk_c", meta$Plate)] <- "booster_1:1:1"
meta$SampleType[meta$Channel == "126" & !grepl("bulk_c", meta$Plate)] <- "booster_bulk"
meta$SampleType[meta$Channel == "127N"] <- "norm"
meta$SampleType[meta$SortedPopulation == "BULK"] <- "sc"
meta$SampleType <- ifelse(meta$SampleType == "empty" & meta$Channel != "127C",
                          "neg control", meta$SampleType)
meta$nbCells <- 0
meta$nbCells[meta$SampleType == "sc"] <- 1
meta$nbCells[meta$SampleType == "norm"] <- 10
meta$nbCells[meta$SampleType == "booster"] <- 200

## Data cleaning
## Some columns are converted to factors
meta$Col <- factor(meta$Col, levels = 1:24)
meta$File.ID <- factor(meta$File.ID, levels = paste0("F", 1:192))
channelOrder <- c("126", "127N", "127C", "128N", "128C", "129N", "129C",
                  "130N", "130C", "131N", "131C", "132N", "132C", "133N",
                  "133C", "134N")
meta$Channel <- factor(meta$Channel, levels = channelOrder)


####---- PSM data ----####

## Load the PSM quantitative data
list.files(path = dataDir, pattern = "PSMs.txt", full.names = TRUE) %>%
    read.table(header = TRUE, sep = "\t") %>%
    ## Rename columns with the TMT label names only
    setNames(gsub("Abundance.", "", colnames(.))) ->
    psms

## Include contaminant info
cont <- read.table(paste0(dataDir, "../contaminants.txt"), header = TRUE)   
psms$isContaminant <- sapply(psms$Master.Protein.Accessions, 
                             function(x) any(strsplit(x, "; ")[[1]] %in% cont$Accession))

## Create the QFeatures object 
schoof2021 <- readSCP(featureData = psms,
                      colData = meta, 
                      batchCol = "File.ID",
                      channelCol = "Channel")


####---- Protein data ----####


## The protein data was loaded at the start of the script as `prot`
## Adapt column names so they match with PSM data
colnames(prot) %>%
    sub(pattern = "Abundances.Grouped.", replacement = "") %>%
    sub(pattern = "[.]", replacement = "") ->
    colnames(prot)

## Remove unnecessary columns
prot <- prot[, !grepl("^Foundin", colnames(prot))]

## Create protein rownames as in SCeptre, that is use the gene symbol,
## use the accession if there is no gene symbol and append the accession
## for duplicated gene symbols
rown <- prot$GeneSymbol
isEmpty <- which(rown == "")
rown[isEmpty] <- prot$Accession[isEmpty]
isDupl <- rown %in% rown[duplicated(rown)]
rown[isDupl] <- paste0(prot$GeneSymbol[isDupl], "_", prot$Accession[isDupl])
prot$rowname <- rown

## Include contaminant info
prot$isContaminant <- prot$Accession %in% cont$Accession

## Convert the protein data to a SCE object
prot <- readSingleCellExperiment(prot, 
                                 ecol = grep("^F\\d", colnames(prot)),
                                 fnames = "rowname")

## Add protein data to the QFeatures
schoof2021 <- addAssay(schoof2021, 
                       prot, 
                       name = "proteins")
## NOTE: could not map psms to proteins because some psms have multiple
## associated proteins (`Master.Protein.Accessions`). This is not yet
## supported in `QFeatures`.


####---- Processed protein data ----####


## We also include the protein data that has been processed by the 
## authors. These data are available from the ZIP file distributed by
## the others and is called `SCeptre/Schoof_et_al/results/bulk/bulk.h5ad`.
## Another way to retrieve this data is to run the notebook
## `SCeptre/Schoof_et_al/results/bulk/bulk.h5ad`. 

## The datasets provided by the authors is loaded from Python using 
## the reticulate interface
library(reticulate)
## To make this work, you need to install the sceptre conda 
## environment as recommended by the author's Github page, see 
## https://github.com/bfurtwa/SCeptre#installation-from-github,
## and activate the environment on the system. Once active, you can 
## call `use_condaenv`.
use_condaenv("sceptre")
dataDir <- normalizePath(dataDir)
py_run_string("
import scanpy
adata = scanpy.read_h5ad(r.dataDir + '/../../results/bulk/bulk.h5ad')
              ")
## The dataset is converted to a SingleCellExperiment object
library(zellkonverter)
logNormProt <- AnnData2SCE(py$adata)

## We match the sample names and proteins
colnames(logNormProt) <- paste0(logNormProt$`File ID`, logNormProt$Channel)

## Keep only the colData columns specific to logNormProt
sel <- c("leiden", "Root Cell", "dpt_pseudotime", "Cluster 5", 
         "Cluster 3", "Cluster 4", "Cluster 2", "Cluster 1")
colData(logNormProt) <- colData(logNormProt)[, sel]

## Insert the data in the `QFeatures` object
schoof2021 <- addAssay(schoof2021, 
                       logNormProt, 
                       name = "logNormProteins")
schoof2021 <- addAssayLink(schoof2021, 
                           from = "proteins",
                           to = "logNormProteins",
                           varFrom = "Accession",
                           varTo = "Accession")


####---- Save dataset ----####


## Save the `QFeatures` object in an Rda file that will be uploaded to
## ExperimentHub
save(schoof2021,
     compress = "xz", 
     compression_level = 9,
     file = "~/PhD/.localdata/scpdata/schoof2021.Rda")

