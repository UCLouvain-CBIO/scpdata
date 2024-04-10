
####---- Petrosius et al, 2023 ---####


## Petrosius, V., Aragon-Fernandez, P., Üresin, N. et al. Exploration of cell
## state heterogeneity using single-cell proteomics through sensitivity-tailored
## data-independent acquisition. Nat Commun 14, 5910 (2023). 
## https://doi.org/10.1038/s41467-023-41602-1

library(SingleCellExperiment)
library(scp)
library(tidyverse)

####---- Load PSM data ----####
## The PSM data downloaded from the https://dataverse.uclouvain.be/dataset.xhtml?persistentId=doi:10.14428/DVN/EMAVLT 
## and 'sample_facs.csv' from the https://zenodo.org/records/8146605
## '20240205_111251_PEPQuant (Normal).tsv' = contains the PSM data.
## 'sample_facs.csv' = contains the cell annotations.

root <- "~/localdata/SCP/petrosiusmESC/20240205_111248_mESC_SNEcombine_m15-m2i/"
ev <- read.delim(paste0(root, "20240205_111251_PEPQuant (Normal).tsv"))
design <- read.delim(paste0(root, "sample_facs.csv"))

####---- Create sample annotation ----####
design %>%
  select(-X) %>%
  distinct() %>%
  add_column(Channel = "PEP.Quantity") %>%
  rename(Set = File.Name, 
         SampleType = Plate) ->
  meta

## Clean quantitative data
ev %>%
  rename(Set = R.FileName, 
         protein = PG.ProteinAccessions) %>%
  ## Create a modified sequence + charge variable
  mutate(peptide = paste0("_", PEP.StrippedSequence, "_.", FG.Charge)) %>%
  filter(Set %in% meta$Set) ->
  evproc

## Create the QFeatures object
petrosius2023_mES <- readSCP(evproc, 
                    meta, 
                    channelCol = "Channel", 
                    batchCol = "Set",
                    removeEmptyCols = TRUE)


####---- Peptide data ----####
## The peptide data downloaded from the https://dataverse.uclouvain.be/dataset.xhtml?persistentId=doi:10.14428/DVN/EMAVLT  
## '20240205_111251_Peptide Quant (Normal).tsv' contains the peptide data.

## Load the peptide level quantification data
pep_data <- read.delim(paste0(root, "20240205_111251_Peptide Quant (Normal).tsv"))

## Clean quantitative data
pep_data %>%
  pivot_wider(names_from = R.FileName, 
              values_from = PG.Quantity, 
              id_cols = c(EG.PrecursorId, PG.ProteinAccessions)) ->
  peps

## Create the SingleCellExperiment object
pep <- readSingleCellExperiment(peps, 
                                ecol = 3:605)

## Name rows with peptide sequence
rownames(pep) <- peps$EG.PrecursorId

## Rename columns so they math with the PSM data
colnames(pep) %>%
  paste0("PEP.Quantity") ->
  colnames(pep)

## Include the peptide data in the QFeatures object
petrosius2023_mES <- addAssay(petrosius2023_mES, pep, name = "peptides")

## Link the PSMs and the peptides
petrosius2023_mES <- addAssayLink(petrosius2023_mES, 
                           from = 1:603, 
                           to = "peptides",
                           varFrom = rep("EG.PrecursorId", 603),
                           varTo = "EG.PrecursorId")


####---- Add the protein data ----####
## The peptide data downloaded from the https://dataverse.uclouvain.be/dataset.xhtml?persistentId=doi:10.14428/DVN/EMAVLT
## '20240205_111251_PGQuant (Normal).tsv' contains the protein data.

prot_data <- read.delim(paste0(root, "20240205_111251_PGQuant (Normal).tsv"))

## Clean quantitative data
prot_data %>% 
  mutate(R.FileName = sub(".*rawfiles/", "", R.Raw.File.Name)) %>%
  mutate(R.FileName = sub(".raw", "", R.FileName)) %>%
  pivot_wider(names_from = R.FileName, 
              values_from = PG.Quantity, 
              id_cols = PG.ProteinAccessions) ->
  prots

## Create the SingleCellExperiment object
pro <- readSingleCellExperiment(prots, 
                                ecol = 2:604)

## Name rows with peptide sequence
rownames(pro) <- prots$PG.ProteinAccessions

## Rename columns so they math with the PSM data
colnames(pro) %>%
  paste0("PEP.Quantity") ->
  colnames(pro)

## Include the peptide data in the QFeatures object
petrosius2023_mES <- addAssay(petrosius2023_mES, pro, name = "proteins")

## Link the PSMs and the peptides
petrosius2023_mES <- addAssayLink(petrosius2023_mES, 
                                from = "peptides", 
                                to = "proteins",
                                varFrom = "PG.ProteinAccessions",
                                varTo = "PG.ProteinAccessions")

## Save data
save(petrosius2023_mES,
     file = file.path(paste0(root, "petrosius2023_mES.Rda")),
     compress = "xz",
     compression_level = 9)

