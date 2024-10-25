
####---- Hu et al, 2023 ---####


## Hu, M., Zhang, Y., Yuan, Y., Ma, W., Zheng, Y., Gu, Q., & Xie, X. S. 2023. 
## “Correlated protein modules revealing functional coordination of interacting 
## proteins are detected by single-cell proteomics.”. The Journal of Physical 
## Chemistry B, https://doi.org/10.1021/acs.jpcb.3c00014 

library(SingleCellExperiment)
library(scp)
library(tidyverse)

root <- "~/localdata/SCP/hu2023/"

####---- Add the protein data ----####

## Data accessible at GitHub repository
## https://github.com/dionezhang/CPM/blob/master/ProteinAbundance.Rdata

#### Load data ####
load(paste0(root, "ProteinAbundance.Rdata"))

Norm %>% 
  mutate(X = rownames(Norm)) %>% 
  readSingleCellExperiment(ecol = 1:69, fnames = "X") ->
  K562

## Protein data for K562 cells
hu2023_K562 <- SingleCellExperiment(K562)

prots <- rownames(hu2023_K562)
rowData(hu2023_K562) <- Description[prots, ,drop = FALSE]
rowData(hu2023_K562)$protein <- prots

colData(hu2023_K562) <- DataFrame(row.names = colnames(Norm), 
                           SampleType = rep("K562", length(colnames(Norm))))

## Save data
save(hu2023_K562,
     file = file.path(paste0(root, "hu2023_K562.Rda")),
     compress = "xz",
     compression_level = 9)
