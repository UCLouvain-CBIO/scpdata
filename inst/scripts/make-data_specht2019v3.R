
####---- Specht et al. 2019 ---####


## Specht, Harrison, Edward Emmott, Toni Koller, and Nikolai Slavov. 2019. 
## “High-Throughput Single-Cell Proteomics Quantifies the Emergence of Macrophage 
## Heterogeneity.” bioRxiv. https://doi.org/10.1101/665307.

## NOTE: the script is version 2 because the data set was updated by the authors

library(SingleCellExperiment)
library(scp)
library(tidyverse)
setwd("../.localdata/SCP/specht2019/v3/")

####---- Load PSM data ----####

## The `raw.RData` files was downloaded from 
## https://drive.google.com/drive/folders/1Zhjik_JFjCQNIVjg63-fooJ4K0HZxWjV
## It contains 3 datasets
## ev = the MaxQuant output file for all batches with the PSM data
## design = the cell annotation
## batch = the batch annotation 
load("raw.RData")

## Clean the sample metadata so that it meets the requirements for
## `scp::readSCP`. The cell annotation and batch annotation are merge into a 
## table
inner_join(x = design %>% 
               pivot_longer(-Set, 
                            names_to = "Channel", 
                            values_to = "SampleAnnotation") %>%
               mutate(SampleType = recode(SampleAnnotation, 
                                          sc_0 = "Blank",
                                          sc_u = "Monocyte",
                                          sc_m0 = "Macrophage",
                                          unused = "Unused",
                                          norm = "Reference",
                                          reference = "Reference",
                                          carrier_mix = "Carrier")), 
           y = batch %>% 
               rename(Set = set) %>%
               mutate_all(as.character),
           by = "Set") %>%
    mutate(Set = sub("^X", "", Set)) ->
    meta

## Clean quantitative data
ev %>%
    ## Variable names should be consistent with metadata
    rename_all(gsub, 
               pattern = "^Reporter[.]intensity[.](\\d*)$", 
               replacement = "RI\\1") %>%
    rename(Set = Raw.file, 
           protein = Leading.razor.protein) %>%
    ## Remove "X" at start of batch 
    mutate(Set = gsub("^X", "", Set),
           ## Create a modified sequence + charge variable
           peptide = paste0(Modified.sequence, Charge)) %>%
    ## keep only single cell runs 
    filter(!grepl("col\\d{2}|arrier|Ref|Master|SQC|blank", Set)) %>%
    ## remove experimental sets concurrent with low mass spec performance
    filter(!grepl("FP9[56]|FP103", Set)) %>%
    ## Make sure all runs are described in design, if not, remove them
    filter(Set %in% meta$Set) ->
    evproc


## Create the QFeatures object
specht2019v3 <- readSCP(evproc, 
                        meta, 
                        channelCol = "Channel", 
                        batchCol = "Set",
                        removeEmptyCols = TRUE)

####---- Get the SCoPE2 annotation table ----####

## The protein and peptide data contain features by single-cells. The 
## columns are named after an arbitrary indexing that must be 
## converted before adding the assay to the QFeatures object to 
## preserve consistency across assays.

## The `Cells.csv` file was downloaded from 
## https://scope2.slavovlab.net/docs/data and contains the cell IDs 
## and annotations. Some processing is preserve consistency with the 
## PSM data.
read.csv("Cells.csv", 
         row.names = 1) %>%
    t %>%
    as.data.frame %>%
    rownames_to_column("id") %>%
    rename(Set = raw.file, 
           SampleAnnotation = celltype, 
           digest = batch_digest, 
           sortday = batch_sort, 
           lcbatch = batch_chromatography) %>%
    mutate(Channel = NA, 
           Set = sub("^X", "", Set),
           SampleType = recode(SampleAnnotation, 
                               sc_0 = "Blank",
                               sc_u = "Monocyte",
                               sc_m0 = "Macrophage",
                               unused = "Unused",
                               norm = "Reference",
                               reference = "Reference",
                               carrier_mix = "Carrier")) %>%
    column_to_rownames("id") %>%
    as("DataFrame") ->
    annot
## The `IDtoChannel.csv` file was generated using the `SCoPE2_analysis.R`
## script from https://github.com/cvanderaa/SCoPE2/tree/master/code.
## The csv contains the data of the `c2q` variable that maps the cell
## id to the RI channel
read.csv("IDtoChannel.csv", )[, -1] %>%
    filter(celltype %in% rownames(annot)) %>%
    mutate(channel = sub("Reporter[.]intensity[.]", "RI", channel)) ->
    idMap
## Add channel info 
annot[idMap$celltype, "Channel"] <- idMap$channel


####---- Add the peptide data ----####


## The `Peptides-raw.csv` file was downloaded from 
## https://scope2.slavovlab.net/docs/data 
pep <- readSingleCellExperiment("Peptides-raw.csv", 
                                ecol = -c(1,2), 
                                fnames = "peptide")
colData(pep) <- annot[colnames(pep), ]
colnames(pep) <- paste0(annot[colnames(pep), "Set"],  annot[colnames(pep), "Channel"])
## Include the peptide data in the QFeatures object
specht2019v3 <- addAssay(specht2019v3, pep, name = "peptides")

## Link the PSMs and the peptides
## First find which PSM assays were included
sel <- sapply(grep("19", names(specht2019v3), value = TRUE), function(name) {
    x <- specht2019v3[[name]]
    ## Does the current PSM data have at least 1 colname in common with pep?
    inColnames <- any(colnames(x) %in% colnames(pep))
    ## Does the current PSM data have at least 1 peptide sequence in common with pep?
    inSequence <- any(rowData(x)$peptide %in% rowData(pep)$peptide)
    return(inColnames && inSequence) ## The PSM assay must fulfill both conditions
})
## Add an AssayLink that bridges the PSM assays and the peptide assay
specht2019v3 <- addAssayLink(specht2019v3, from = which(sel), to = "peptides", 
                             varFrom = rep("peptide", sum(sel)), varTo = "peptide")

####---- Add the protein data ----####

## The `Proteins-processed.csv` and `Cells.csv` file was downloaded from 
## https://scope2.slavovlab.net/docs/data
read.csv("Proteins-processed.csv") %>%
    select(-X) %>%
    readSingleCellExperiment(ecol = 1:1490, fnames = "protein") ->
    prot
colData(prot) <- annot[colnames(prot), ]
colnames(prot) <- paste0(annot[colnames(prot), "Set"],  
                         annot[colnames(prot), "Channel"])
specht2019v3 <- addAssay(specht2019v3, prot, name = "proteins")
specht2019v3 <- addAssayLink(specht2019v3, from = "peptides", to = "proteins", 
                             varFrom = "protein", varTo = "protein")

## Save data
save(specht2019v3,
     file = file.path("~/PhD/.localdata/scpdata/specht2019v3.Rda"),
     compress = "xz", 
     compression_level = 9)

