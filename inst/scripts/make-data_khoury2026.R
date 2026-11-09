## Requred packages of this workflow
library("scp")
library("tidyverse")

## All data files were retrieved from this Google Drive repository:
## https://zenodo.org/records/22649483
dataDir <- "~/localdata/SCP/khoury2026/"

## Sample annotations
metaData <- read.csv(paste0(dataDir, "PBMC_covariation_protein_cell_metadata.csv"))

metaData <- metaData |>
  mutate(runCol = lc_run_name,
         quantCols = mtraq_channel)

# PSM Data
reportFiles <- list.files(paste0(dataDir, "Batch_IDs_and_abundances/"),
                          pattern = "IdentificationsList_", full.names = FALSE)

## Load and rowbind all report files
psms <- do.call(rbind,
                lapply(reportFiles, function(x) {
                  df <- read.csv(x)
                  df$dataset_batch <- sub(".*_(.*?)\\.csv$", "\\1", x)
                  df
                }
                ))

psms$cell_id <- paste0(psms$File.Name, "_", psms$dataset_batch)

psmsWide <- psms |>
  left_join(select(metaData, cell_id, lc_run_name, mtraq_channel),
            by = "cell_id") |>
  select(-File.Name, -cell_id) |>
  pivot_wider(names_from = mtraq_channel, values_from = Ms1.Area) |>
  as.data.frame()

## Format to a QFeatures object
khoury2026 <- readSCP(assayData = psmsWide,
                      colData = DataFrame(metaData),
                      runCol = "lc_run_name",
                      removeEmptyCols = TRUE)

## Map cell IDs to QFeatures sample names
cd <- colData(khoury2026)
sampleNames <- setNames(rownames(cd), cd$cell_id)

# Peptide Data
## Gene to Protein Map
P2G <- read.csv(paste0(dataDir, "PepProtMap_integrated.csv"))
protMap <- distinct(P2G, Protein.Group, Genes)

pep1 <- read.csv(paste0(dataDir, "Peptide_x_Cell_Matrix_Normalized_Centered.csv"))
pep2 <- read.csv(paste0(dataDir, "Peptide_x_Cell_Matrix_Column_Normalized.csv"))

## Create the SingleCellExperiment object
pep1 <- readSingleCellExperiment(pep1, quantCols = 2:ncol(pep1))
colnames(pep1) <- unname(sampleNames[colnames(pep1)])
rownames(pep1) <- rowData(pep1)$Stripped.Sequence
rowData(pep1) <- cbind(rowData(pep1),
                       P2G[match(rownames(pep1), P2G$Stripped.Sequence), -1])

pep2 <- readSingleCellExperiment(pep2, quantCols = 2:ncol(pep2))
colnames(pep2) <- unname(sampleNames[colnames(pep2)])
rownames(pep2) <- rowData(pep2)$Stripped.Sequence
rowData(pep2) <- cbind(rowData(pep2),
                       P2G[match(rownames(pep2), P2G$Stripped.Sequence), -1])

from <- names(khoury2026)
khoury2026 <- addAssay(khoury2026, pep1, "Peptide_Normalized_Centered")
khoury2026 <- addAssayLink(khoury2026,
                           from = from,
                           to = "Peptide_Normalized_Centered",
                           varFrom = rep("Stripped.Sequence", length(from)),
                           varTo = "Stripped.Sequence")

khoury2026 <- addAssay(khoury2026, pep2, "Peptide_Column_Normalized")
khoury2026 <- addAssayLink(khoury2026,
                           from = from,
                           to = "Peptide_Column_Normalized",
                           varFrom = rep("Stripped.Sequence", length(from)),
                           varTo = "Stripped.Sequence")

## Protein Data
prot1 <- read.csv(paste0(dataDir, "Protein_x_Cell_Matrix_Imputed_Column_Normalized.csv"))
prot2 <- read.csv(paste0(dataDir, "Protein_x_Cell_Matrix_Imputed_Normalized_Centered.csv"))
prot3 <- read.csv(paste0(dataDir, "Protein_x_Cell_Matrix_Normalized_Centered.csv"))
prot4 <- read.csv(paste0(dataDir, "Protein_x_Cell_Matrix_Column_Normalized.csv"))

## Create the SingleCellExperiment objects
prot1 <- readSingleCellExperiment(prot1, quantCols = 2:ncol(prot1))
colnames(prot1) <- unname(sampleNames[colnames(prot1)])
rownames(prot1) <- rowData(prot1)$Genes
rowData(prot1)$Protein.Group <- protMap$Protein.Group[match(rownames(prot1),
                                                            protMap$Genes)]

prot2 <- readSingleCellExperiment(prot2, quantCols = 2:ncol(prot2))
colnames(prot2) <- unname(sampleNames[colnames(prot2)])
rownames(prot2) <- rowData(prot2)$Genes
rowData(prot2)$Protein.Group <- protMap$Protein.Group[match(rownames(prot2),
                                                            protMap$Genes)]

prot3 <- readSingleCellExperiment(prot3, quantCols = 2:ncol(prot3))
colnames(prot3) <- unname(sampleNames[colnames(prot3)])
rownames(prot3) <- rowData(prot3)$Genes
rowData(prot3)$Protein.Group <- protMap$Protein.Group[match(rownames(prot3),
                                                            protMap$Genes)]

prot4 <- readSingleCellExperiment(prot4, quantCols = 2:ncol(prot4))
colnames(prot4) <- unname(sampleNames[colnames(prot4)])
rownames(prot4) <- rowData(prot4)$Genes
rowData(prot4)$Protein.Group <- protMap$Protein.Group[match(rownames(prot4),
                                                            protMap$Genes)]

## Add protein assays to QFeatures and link
khoury2026 <- addAssay(khoury2026, prot1, "Protein_Imputed_Column_Normalized")
khoury2026 <- addAssayLink(khoury2026,
                           from = "Peptide_Column_Normalized",
                           to = "Protein_Imputed_Column_Normalized",
                           varFrom = "Genes",
                           varTo = "Genes")

khoury2026 <- addAssay(khoury2026, prot2, "Protein_Imputed_Normalized_Centered")
khoury2026 <- addAssayLink(khoury2026,
                           from = "Peptide_Normalized_Centered",
                           to = "Protein_Imputed_Normalized_Centered",
                           varFrom = "Genes",
                           varTo = "Genes")

khoury2026 <- addAssay(khoury2026, prot3, "Protein_Normalized_Centered")
khoury2026 <- addAssayLink(khoury2026,
                           from = "Peptide_Normalized_Centered",
                           to = "Protein_Normalized_Centered",
                           varFrom = "Genes",
                           varTo = "Genes")

khoury2026 <- addAssay(khoury2026, prot4, "Protein_Column_Normalized")
khoury2026 <- addAssayLink(khoury2026,
                           from = "Peptide_Column_Normalized",
                           to = "Protein_Column_Normalized",
                           varFrom = "Genes",
                           varTo = "Genes")

## Save data as rda file
save(khoury2026,
     file = paste0(dataDir, "khoury2026.rda"),
     compress = "xz",
     compression_level = 9)

