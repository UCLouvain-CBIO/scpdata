
library("scp")
library("tidyverse")

## library("devtools")
## devtools::install_github("https://github.com/SlavovLab/QuantQC")
library("QuantQC") ## Required to access QuantQC data slots @

dir <- "/storage/research/dduv/cbio-lg/user/eayar/localdata/SCP/leduc2025/dat/"

# Processed peptide level data
load(paste0(dir, "03_QuantQC_objects/r1_5day_male.RData"))
load(paste0(dir, "03_QuantQC_objects/r2_5day_female.RData"))
load(paste0(dir, "03_QuantQC_objects/r3_10day_male.RData"))
load(paste0(dir, "03_QuantQC_objects/r4_10day_female.RData"))

## Extract peptide to protein map from data
make_pep2prot <- function(raw_data) {
  raw_data |>
    distinct(seqcharge, Stripped.Sequence, Protein.Group, Genes) |>
    as.data.frame()
}

map1 <- make_pep2prot(r1_5day_male@raw_data)
map2 <- make_pep2prot(r2_5day_female@raw_data)
map3 <- make_pep2prot(r3_10day_male@raw_data)
map4 <- make_pep2prot(r4_10day_female@raw_data)

## Combine row data with quant data into single cell experiment
make_pep_se <- function(mat, map) {
  rd <- map[match(rownames(mat), map$seqcharge), ]
  rownames(rd) <- rownames(mat)

  SummarizedExperiment(
    assays = list(quantification = as.matrix(mat)),
    rowData = DataFrame(rd)
  )
}

r1_pep <- r1_5day_male@matricies@peptide
r2_pep <- r2_5day_female@matricies@peptide
r3_pep <- r3_10day_male@matricies@peptide
r4_pep <- r4_10day_female@matricies@peptide

## Rename samples so they match to protein data
colnames(r1_pep) <- paste0(colnames(r1_pep), "_prep1")
colnames(r2_pep) <- paste0(colnames(r2_pep), "_prep2")
colnames(r3_pep) <- paste0(colnames(r3_pep), "_prep3")
colnames(r4_pep) <- paste0(colnames(r4_pep), "_prep4")

r1_pep_se <- make_pep_se(r1_pep, map1)
r2_pep_se <- make_pep_se(r2_pep, map2)
r3_pep_se <- make_pep_se(r3_pep, map3)
r4_pep_se <- make_pep_se(r4_pep, map4)

# Raw peptide level data
r1_raw <- read.csv(paste0(dir, "02_raw_reptide_X_singleCell/r1_peptide.csv"), row.names = 1)
r2_raw <- read.csv(paste0(dir, "02_raw_reptide_X_singleCell/r2_peptide.csv"), row.names = 1)
r3_raw <- read.csv(paste0(dir, "02_raw_reptide_X_singleCell/r3_peptide.csv"), row.names = 1)
r4_raw <- read.csv(paste0(dir, "02_raw_reptide_X_singleCell/r4_peptide.csv"), row.names = 1)

## Rename samples so they match to protein data
colnames(r1_raw) <- paste0(colnames(r1_raw), "_prep1")
colnames(r2_raw) <- paste0(colnames(r2_raw), "_prep2")
colnames(r3_raw) <- paste0(colnames(r3_raw), "_prep3")
colnames(r4_raw) <- paste0(colnames(r4_raw), "_prep4")

r1_raw_se <- make_pep_se(r1_raw, map1)
r2_raw_se <- make_pep_se(r2_raw, map2)
r3_raw_se <- make_pep_se(r3_raw, map3)
r4_raw_se <- make_pep_se(r4_raw, map4)

## Generate QFeatures object with peptide level data
leduc2025 <- QFeatures(list(
  pep_raw_prep1 = r1_raw_se,
  pep_raw_prep2 = r2_raw_se,
  pep_raw_prep3 = r3_raw_se,
  pep_raw_prep4 = r4_raw_se,
  pep_processed_prep1 = r1_pep_se,
  pep_processed_prep2 = r2_pep_se,
  pep_processed_prep3 = r3_pep_se,
  pep_processed_prep4 = r4_pep_se
))

# Protein Level shared data
## Absolute abundances
absolute <- read.csv(paste0(dir, "04_Gene_X_SingleCell_and_annotations/sc_protein_absolute.csv"),
                     row.names = 1)
prot_absolute_se <- SummarizedExperiment(assay = as.matrix(absolute),
                                         rowData = DataFrame(Protein = rownames(absolute)))

## Relative abundances
relative <- read.csv(paste0(dir, "04_Gene_X_SingleCell_and_annotations/sc_protein_relative.csv"),
                     row.names = 1)
prot_relative_se <- SummarizedExperiment(assay = as.matrix(relative),
                                         rowData = DataFrame(Protein = rownames(relative)))

## Imputed and comBat corrected data used in the study
imputed <- readRDS(paste0(dir, "protein_Data.rds"))
prot_imputed_se <- SummarizedExperiment(assay = as.matrix(imputed),
                                         rowData = DataFrame(Protein = rownames(imputed)))

## colData object
annot <- read.csv(paste0(dir, "04_Gene_X_SingleCell_and_annotations/sc_protein_annotations.csv"))
rownames(annot) <- annot$ID

## Add protein level assays
leduc2025 <- addAssay(leduc2025, prot_relative_se, name = "proteins_relative")
leduc2025 <- addAssay(leduc2025, prot_absolute_se, name = "proteins_absolute")
leduc2025 <- addAssay(leduc2025, prot_imputed_se, name = "proteins_imputed")

## Add sample meta data
colData(leduc2025) <- DataFrame(annot)

## Link corresponding peptide assays to protein level
leduc2025 <- addAssayLink(leduc2025,
                          from = c("pep_raw_prep1", "pep_raw_prep2",
                                   "pep_raw_prep3", "pep_raw_prep4"),
                          to = "proteins_absolute",
                          varFrom = rep("Protein.Group", 4),
                          varTo = "Protein")
leduc2025 <- addAssayLink(leduc2025,
                          from = c("pep_processed_prep1", "pep_processed_prep2",
                                   "pep_processed_prep3", "pep_processed_prep4"),
                          to = "proteins_relative",
                          varFrom = rep("Protein.Group", 4),
                          varTo = "Protein")

## Save data as Rda file
save(leduc2025,
     file = paste0(dir, "leduc2025.rda"),
     compress = "xz",
     compression_level = 9)

