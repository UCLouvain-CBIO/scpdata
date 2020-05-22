
####---- Specht et al. 2019 ---####


## Specht, Harrison, Edward Emmott, Toni Koller, and Nikolai Slavov. 2019. 
## “High-Throughput Single-Cell Proteomics Quantifies the Emergence of Macrophage 
## Heterogeneity.” bioRxiv. https://doi.org/10.1101/665307.

## NOTE: the script is version 2 because the data set was updated by the authors

library(Features)
library(SingleCellExperiment)
library(tidyverse)
setwd("./inst/scripts")

## The data was downloaded from https://scope2.slavovlab.net/docs/data to 
## scpdata/inst/extdata/specht2019-v2

## Sample metadata 
## https://drive.google.com/file/d/1c0nUlgQN1CjeWMZGjC_bTmcAe-e9rfaT/view
read.csv("../extdata/specht2019-v2/Cells.csv", 
         header = TRUE, row.names = 1) %>%
  t %>% 
  data.frame -> meta

## Peptide data 
## https://drive.google.com/open?id=1c9fwYI4qf9LzaQHf0wXoYAr2fuqNhNll
pep <- read.csv("../extdata/specht2019-v2/Peptides-raw.csv")
rownames(pep) <- pep$peptide
pep <- SummarizedExperiment(assays = SimpleList(exprs = pep %>%
                                                  select(starts_with("i")) %>%
                                                  as.matrix),
                            rowData = pep %>%
                              select(starts_with("p")) %>%
                              DataFrame)
meta[colnames(pep), ] %>%
  DataFrame -> colData(pep)

## Protein data
## https://drive.google.com/open?id=1c5Z3b_2gOwDyHCLm9ycY3hXckY1GDd5L
prot <- read.csv("../extdata/specht2019-v2/Proteins-processed.csv", 
                 row.names = 1) 
prot <- SummarizedExperiment(assays = SimpleList(exprs = prot %>%
                                                   as.matrix),
                             rowData = DataFrame(protein = rownames(prot)))
meta[colnames(prot), ] %>%
  DataFrame -> colData(prot)

## Create Features object 
Features(ExperimentList(peptide = pep,
                        protein = prot)) %>%
  ## Link peptide assay and protein assay
  addAssayLink(from = "peptide", to = "protein", varFrom = "protein", 
               varTo = "protein") -> specht2019v2_fts

## Test the Features object 
specht2019v2_fts["Q86U42",,] %>%
  longFormat %>%
  data.frame %>%
  mutate(isNa = is.na(value)) %>%
  replace_na(list(value = -2)) %>%
  ggplot(aes(y = value, x = rowname)) +
  geom_violin(aes(col = assay)) +
  geom_jitter(alpha = 0.5) +
  facet_grid(~ assay, scales = "free_x", space = "free_x") +
  theme(legend.position = "none")


## Save data
save(specht2019v2_fts, file = file.path("../../data/specht2019v2_fts.rda"),
     compress = "xz", compression_level = 9)

