
####---- Khan et al, 2023 ---####


## Saad Khan, Rachel Conover, Anand R. Asthagiri, Nikolai Slavov. 2023. 
## “Dynamics of single-cell protein covariation during epithelial–mesenchymal 
## transition.” bioRxiv. https://doi.org/10.1101/2023.12.21.572913.

library(SingleCellExperiment)
library(scp)
library(tidyverse)

####---- Load PSM data ----####

## 3 Seperate files downloaded from the drive link:
## https://drive.google.com/drive/folders/1zCsRKWNQuAz5msxx0DfjDrIe6pUjqQmj
## 'ev_updated_NS.DIA.txt' = the MaxQuant output file for all batches with the PSM data
## 'annotation.csv' = the cell annotation (remainder operation needs to be used //16)
## 'batch.csv' = the batch annotation 

root <- "~/localdata/SCP/khan2023/002-singleCellDataGeneration/"
ev <- read.delim(paste0(root, "ev_updated_NS.DIA.txt"))
design <- read.csv(paste0(root, "annotation.csv"))
batch <- read.csv(paste0(root, "batch.csv"))

## Clean the sample metadata so that it meets the requirements for
## `scp::readSCP`. The cell annotation and batch annotation are merged into a 
## table
inner_join(x = design %>% 
             pivot_longer(-Set, 
                          names_to = "Channel", 
                          values_to = "SampleAnnotation") %>%
             mutate(SampleType = recode(SampleAnnotation,
                                        neg = "Negative",
                                        d0 = "Day_0",
                                        d3 = "Day_3",
                                        d9 = "Day_9",
                                        unused = "Unused",
                                        reference = "Reference",
                                        carrier = "Carrier")), 
           y = batch %>% 
             rename(Set = set) %>%
             mutate_all(as.character),
           by = "Set") -> meta

## Clean quantitative data
ev %>%
  ## Variable names should be consistent with metadata
  rename_all(gsub, 
             pattern = "^Reporter[.]intensity[.](\\d*)$", 
             replacement = "RI\\1") %>%
  rename(Set = Raw.file, 
         protein = Leading.razor.protein) %>%
  ## Create a modified sequence + charge variable
  mutate(peptide = paste0(Modified.sequence, Charge)) %>%
  ## keep only single cell runs 
  filter(!grepl("col\\d{2}|arrier|Ref|Master|SQC|blank", Set)) %>%
  ## remove experimental sets concurrent with low mass spec performance
  filter(!grepl("FP9[56]|FP103", Set)) %>%
  ## Make sure all runs are described in design, if not, remove them
  filter(Set %in% meta$Set) ->
  evproc

## Create the QFeatures object
khan2023 <- readSCP(evproc, 
                        meta, 
                        channelCol = "Channel", 
                        batchCol = "Set",
                        removeEmptyCols = TRUE)

####---- Get the annotation table ----####

## The protein data contain features by single-cells. The 
## columns are named after an arbitrary indexing that must be 
## converted before adding the assay to the QFeatures object to 
## preserve consistency across assays.

## Annotation file was downloaded from 
## https://drive.google.com/drive/folders/1zCsRKWNQuAz5msxx0DfjDrIe6pUjqQmj
## 'annotation.csv': Cell annotations

## The `IDtoChannel.csv` file was generated using the `EMTTGFB_singleCellProcessing.R`
## script from https://github.com/SlavovLab/EMT_TGFB_2023/tree/main adjusted to
## obtain mapping from the cell id (i...) to the (Set + RI channel)

read.csv(paste0(root, "annotation.csv")) %>%
  pivot_longer(-Set, 
               names_to = "Channel", 
               values_to = "SampleAnnotation") %>%
  mutate(SampleType = recode(SampleAnnotation, 
                             neg = "Negative",
                             d0 = "Day_0",
                             d3 = "Day_3",
                             d9 = "Day_9",
                             unused = "Unused",
                             reference = "Reference",
                             carrier = "Carrier")) %>%
  add_column(lcbatch = "A", sortday = "B", digest = "C") %>%
  unite(rowname, Set, Channel, sep = "", remove = FALSE) %>%
  column_to_rownames() ->
  annot

idMap <- read.csv(paste0(root, "cellIDToChannel.csv"), row.names = 1)

####---- Add the peptide data ----####

## Peptide quantity matrix downloaded from:  
## https://drive.google.com/drive/folders/1zCsRKWNQuAz5msxx0DfjDrIe6pUjqQmj

peps <- read.delim(paste0(root, "EpiToMesen.TGFB.nPoP_trial1_pepByCellMatrix_NSThreshDART_medIntCrNorm.txt"))
peps %>%
  rename(peptide = pep) %>%
  readSingleCellExperiment(ecol = 1:421, fnames = "peptide") ->
  peptides

colnames(peptides) <- idMap$Channel[match(colnames(peptides), idMap$cellID)]
colData(peptides) <- DataFrame(annot[colnames(peptides), ])

khan2023 <- addAssay(khan2023, peptides, name = "peptides")

## Include rowData to peptides assay
rowData(khan2023[["peptides"]]) <- DataFrame(peptide = peps$pep, 
                                             protein = peps$prot)

## First find which PSM assays were included
sel <- sapply(grep("eSK", names(khan2023), value = TRUE), 
              function(name) {
                x <- khan2023[[name]]
                ## Does the current PSM data have at least 1 colname in common with pep?
                inColnames <- any(colnames(x) %in% colnames(peptides))
                ## Does the current PSM data have at least 1 peptide sequence in common with pep?
                inSequence <- any(rowData(x)$peptide %in% rowData(peptides)$peptide)
                return(inColnames && inSequence) ## The PSM assay must fulfill both conditions
})

## Add an AssayLink that bridges the PSM assays and the peptide assay
khan2023 <- addAssayLink(khan2023, from = which(sel), to = "peptides", 
                             varFrom = rep("peptide", sum(sel)), varTo = "peptide")

####---- Add the protein data ----####

## Imputed and un-imputed protein quantity matrices downloaded from:  
## https://drive.google.com/drive/folders/1zCsRKWNQuAz5msxx0DfjDrIe6pUjqQmj

## Add imputed protein data
read.csv(paste0(root, "EpiToMesen.TGFB.nPoP_trial1_ProtByCellMatrix_NSThreshDART_medIntCrNorm_imputedNotBC.csv")) %>%
  rename(protein = X) %>%
  readSingleCellExperiment(ecol = 2:422, fnames = "protein") ->
  prot_imp

colnames(prot_imp) <- idMap$Channel[match(colnames(prot_imp), idMap$cellID)]
colData(prot_imp) <- DataFrame(annot[colnames(prot_imp), ])

khan2023 <- addAssay(khan2023, prot_imp, name = "proteins_imputed")

## Add an AssayLink that bridges the PSM assays and the peptide assay
khan2023 <- addAssayLink(khan2023, from = "peptides", to = "proteins_imputed", 
                         varFrom = "protein", varTo = "protein")

## Add un-imputed protein data
read.csv(paste0(root, "EpiToMesen.TGFB.nPoP_trial1_ProtByCellMatrix_NSThreshDART_medIntCrNorm_unimputed.csv")) %>%
  rename(protein = X) %>%
  readSingleCellExperiment(ecol = 2:422, fnames = "protein") ->
  prot_unimp

colnames(prot_unimp) <- idMap$Channel[match(colnames(prot_unimp), idMap$cellID)]
colData(prot_unimp) <- DataFrame(annot[colnames(prot_unimp), ])

khan2023 <- addAssay(khan2023, prot_unimp, name = "proteins_unimputed")

## Add an AssayLink that bridges the PSM assays and the peptide assay
khan2023 <- addAssayLink(khan2023, from = "peptides", to = "proteins_unimputed", 
                         varFrom = "protein", varTo = "protein")

## Save data
save(khan2023,
     file = file.path(paste0(root, "khan2023.Rda")),
     compress = "xz", 
     compression_level = 9)

