
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

## Data shared by the author, and accessible at 
## https://biopic-my.sharepoint.cn/:x:/g/personal/humo_biopic_pku_edu_cn/EfX4CHedVopLuSx2OJNj6LABdESGNdKz4Eh8Zawvd-fNNQ?rtime=7Xzb4B303Eg

#### Load Data ####
oocyte <- read.csv(paste0(root, "DataMatrix-oocyte-20240614.csv"))
oocyte %>% 
  rename(protein = X) %>% 
  readSingleCellExperiment(ecol = 2:138, fnames = "protein") ->
  oocyte

## Protein data for oocytes
hu2023_oocyte <- SingleCellExperiment(oocyte)

colData(hu2023_oocyte) <- DataFrame(row.names = colnames(hu2023_oocyte), 
                           SampleType = rep("oocyte", length(colnames(oocyte))))

## Save data
save(hu2023_oocyte,
     file = file.path(paste0(root, "hu2023_oocyte.Rda")),
     compress = "xz",
     compression_level = 9)

