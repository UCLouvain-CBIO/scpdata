
####---- Bubis et al. 2025 ---####

libPath <- "/storage/research/dduv/cbio-lg/user/samgregoire/2022-phd-samuel-gregoire/libs/"

library("SingleCellExperiment")
library("QFeatures", lib.loc = libPath)
library("scp", lib.loc = libPath)
library("dplyr")
library("tidyr")
library("purrr")
library("stringr")

# All files were downloaded from the proteomeXchange plateform
# ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2025/06/PXD064564

dataDir <-
  "/storage/research/dduv/cbio-lg/cluster/DataSets/SCPCBIO/bubis2025_astral/"

fileNames <- sub("\\.tsv", "", dir(dataDir))

####---- PSM data ----####

## Read data
PSMFileNames <- fileNames[grepl("Peptides|Precursors", fileNames)]
PSMList <- vector("list", length(PSMFileNames))
PSMList <- map(dir(dataDir, full.names = TRUE)[grepl("Peptides|Precursors", fileNames)],
               read.delim)

names(PSMList) <- PSMFileNames

## Keep only relevant cols
PSMKeepCols <-
  PSMList |>
  map(colnames) |>
  map(grepl, pattern = "EG.PrecursorId|PEP.MS1Quantity")

PSMKeepCols <- map2(map(PSMList, colnames), PSMKeepCols, `[`)
PSMList <- map2(PSMList, PSMKeepCols, select)

## Remove characters from quant cols
for (i in seq_along(PSMList)) {
  PSMList[[i]][PSMList[[i]] == "Filtered"] <- NaN
}

## Convert all quant cols to numeric
for (i in seq_along(PSMList)) {
  PSMList[[i]][, -(which(colnames(PSMList[[i]]) == "EG.PrecursorId"))] <-
    as.data.frame(
      sapply(PSMList[[i]][, -(which(colnames(PSMList[[i]]) == "EG.PrecursorId"))], as.numeric))
}

## Manual annotation of the search strategy
## Based on the README.xslx file

for (i in grep("_MBR|all_ME", names(PSMList))) {
  PSMList[[i]]$searchStrategy <- "directDIA+"
}

for (i in grep("_noMBR|method_eval", names(PSMList))) {
  PSMList[[i]]$searchStrategy <- "method evaluation"
}

for (i in grep("MM_Astral_paper_A549_Figure7_MBR_perGroup_directDIA_strictFDR|JB_250pgK562_50SPD_MBR_OS|JB_250pgHeLa_50SPD_SO_MBR", names(PSMList))) {
  PSMList[[i]]$searchStrategy <- "directDIA+ OS"
}

for (i in grep("JB_HeLa250pg_240k_20Th_40ms_50SPD_L10ng|MM_Astral_paper_K562_50SPD_library_search", names(PSMList))) {
  PSMList[[i]]$searchStrategy <- "library search 10ng"
}

PSMList[[grep("MM_AstralPaper_LibrarySearch_20x", names(PSMList))]]$searchStrategy <-
  "library search 20SC"

PSMList[[grep("MM_AstralPaper_LibrarySearch_40x", names(PSMList))]]$searchStrategy <-
  "library search 40SC"

for (i in grep("L100SC", names(PSMList))) {
  PSMList[[i]]$searchStrategy <- "library search 100SC"
}

## Pivot in long format

PSMList <-
  map(.x = PSMList, .f = pivot_longer,
      names_to = "file", values_to = "PEP.MS1Quantity",
      -c(EG.PrecursorId, searchStrategy))

## Bind list into one dataframe

PSM <- do.call("rbind", PSMList)
rm(PSMList)

## Remove PSMs without quantification

PSM <- PSM[!is.na(PSM$PEP.MS1Quantity), ]

## Clean cols

PSM$Stripped.Sequence <- gsub("_", "", PSM$EG.PrecursorId)
PSM$Stripped.Sequence <- sub("\\..", "", PSM$Stripped.Sequence)
PSM$Stripped.Sequence <- gsub("\\[.*\\]", "", PSM$Stripped.Sequence)

PSM$file <- sub("X\\.\\d+\\.\\.", "", PSM$file)
PSM$file <- sub("\\.raw\\..*", ".raw", PSM$file)

PSM$key <- paste(PSM$file, gsub(" ", "_", PSM$searchStrategy), sep = "_")

####---- Create sample annotation ----####

sa <- data.frame(runCol = unique(PSM$key),
                 cellType = rep(NA, length(unique(PSM$key))),
                 fig = rep(NA, length(unique(PSM$key))))

## File
sa$file <- sub("\\.raw_.*", ".raw", sa$runCol)

## Search strategy
sa$searchStrategy <- sub(".*\\.raw_", "", sa$runCol)

## Cell type
sa$cellType[grepl("HeLa", sa$runCol)] <- "HeLa"
sa$cellType[grepl("A549", sa$runCol)] <- "A549"
sa$cellType[grepl("blank", sa$runCol)] <- "blank"
sa$cellType[grepl("K562", sa$runCol)] <- "K562"
sa$cellType[grepl("H460", sa$runCol)] <- "H460"
sa$cellType[grepl("_SC_TE", sa$runCol)] <- "TE"
sa$cellType[grepl("_SC_hPSC", sa$runCol)] <- "hPSC"

## Throughput
sa$throughput <- str_extract(sa$runCol, "\\d+SPD")

## CV
sa$FAIMS_CV <- "-48V"
sa$FAIMS_CV[grep("FAIMSCV.38", sa$runCol)] <- "-38V"
sa$FAIMS_CV[grep("FAIMSCV.58", sa$runCol)] <- "-58V"
sa$FAIMS_CV[grep("FAIMSCV.68", sa$runCol)] <- "-68V"
sa$FAIMS_CV[grep("FAIMSCV.78", sa$runCol)] <- "-78V"
sa$FAIMS_CV[grep("FAIMSCV.88", sa$runCol)] <- "-88V"
sa$FAIMS_CV[grep("noFAIMS", sa$runCol)] <- "no FAIMS"

## Figure
sa$fig <- NA
sa$fig[grep("DDM|Nov_col", sa$runCol)] <- "fig 1"
sa$fig[grep("_CV\\..8|noFAIMS", sa$runCol)] <- "fig 2"
sa$fig[grep("A549|H460|50SPD_blank|blank_50SPD", sa$file)] <- "fig 4"
sa$fig[grep("gas3p8_blank|gas3p8_SC_TE|gas3p8_SC_hPSC", sa$file)] <- "fig 5"

####---- Build scp object ----####

bubis2025 <- readSCP(
  assayData = PSM,
  colData = sa,
  runCol = "key",
  quantCols = "PEP.MS1Quantity")
rm(PSM)

####---- Generate peptide data ----####

## Aggregate
bubis2025 <-
  aggregateFeatures(bubis2025,
                    i = 1:length(bubis2025),
                    fcol = "Stripped.Sequence",
                    fun = colMedians,
                    name = paste0("peptide_", names(bubis2025)))

## Join assays
bubis2025 <- joinAssays(bubis2025,
                  i = grep("peptide_", names(bubis2025))[bubis2025$fig == "fig 1"],
                  name = "peptides_fig1")
colData(bubis2025[["peptides_fig1"]]) <- colData(bubis2025)[bubis2025$fig == "fig 1", ]

bubis2025 <- joinAssays(bubis2025,
                  i = grep("peptide_", names(bubis2025))[bubis2025$fig == "fig 2"],
                  name = "peptides_fig2")
colData(bubis2025[["peptides_fig2"]]) <- colData(bubis2025)[bubis2025$fig == "fig 2", ]

bubis2025 <- joinAssays(bubis2025,
                  i = grep("peptide_", names(bubis2025))[bubis2025$fig == "fig 4"],
                  name = "peptides_fig4")
colData(bubis2025[["peptides_fig4"]]) <- colData(bubis2025)[bubis2025$fig == "fig 4", ]

bubis2025 <- joinAssays(bubis2025,
                  i = grep("peptide_", names(bubis2025))[bubis2025$fig == "fig 5"],
                  name = "peptides_fig5")
colData(bubis2025[["peptides_fig5"]]) <- colData(bubis2025)[bubis2025$fig == "fig 5", ]

## Remove individual peptide sets

bubis2025 <- bubis2025[, , !grepl("peptide_", names(bubis2025))]

####---- Retrieve Protein data ----####

## Read files
protFileNames <- fileNames[grepl("Protein", fileNames)]
protList <- vector("list", length(protFileNames))
protList <- map(dir(dataDir, full.names = TRUE)[grepl("Protein", fileNames)],
                read.delim)

names(protList) <- protFileNames

## Keep only relevant cols

protKeepCols <-
  protList |>
  map(colnames) |>
  map(grepl, pattern = "PG.Quantity|PG.ProteinGroups")

protKeepCols <- map2(map(protList, colnames), protKeepCols, `[`)

protList <- map2(protList, protKeepCols, select)

## Convert all quant cols to numeric

for (i in seq_along(protList)) {
  protList[[i]][, -(which(colnames(protList[[i]]) == "PG.ProteinGroups"))] <-
    as.data.frame(
      sapply(protList[[i]][, -(which(colnames(protList[[i]]) == "PG.ProteinGroups"))], as.numeric))
}

## Manual annotation of the search strategy
## Based on the README.xslx file

for (i in grep("_MBR|all_ME", names(protList))) {
  protList[[i]]$searchStrategy <- "directDIA+"
}

for (i in grep("_noMBR|method_eval", names(protList))) {
  protList[[i]]$searchStrategy <- "method evaluation"
}

for (i in grep("MM_Astral_paper_A549_Figure7_MBR_perGroup_directDIA_strictFDR|JB_250pgK562_50SPD_MBR_OS|JB_250pgHeLa_50SPD_SO_MBR", names(protList))) {
  protList[[i]]$searchStrategy <- "directDIA+ OS"
}

for (i in grep("JB_HeLa250pg_240k_20Th_40ms_50SPD_L10ng|MM_Astral_paper_K562_50SPD_library_search", names(protList))) {
  protList[[i]]$searchStrategy <- "library search 10ng"
}

protList[[grep("MM_AstralPaper_LibrarySearch_20x", names(protList))]]$searchStrategy <-
  "library search 20SC"

protList[[grep("MM_AstralPaper_LibrarySearch_40x", names(protList))]]$searchStrategy <-
  "library search 40SC"

for (i in grep("L100SC", names(protList))) {
  protList[[i]]$searchStrategy <- "library search 100SC"
}

## Pivot into long format

protList <-
  map(.x = protList, .f = pivot_longer,
      names_to = "file", values_to = "PG.Quantity",
      -c(PG.ProteinGroups, searchStrategy))

## Bind into one dataframe

prot <- do.call("rbind", protList)
rm(protList)

## Remove PSMs without quantification

prot <- prot[!is.na(prot$PG.Quantity), ]

## Clean cols

prot$file <- sub("X\\.\\d+\\.\\.", "", prot$file)
prot$file <- sub("\\.raw\\..*", ".raw", prot$file)

prot$key <- paste(prot$file, gsub(" ", "_", prot$searchStrategy), sep = "_")

####---- Build scp object for proteins ----####

## readSCP
bubis2025_prot <- readSCP(assayData = prot,
                          colData = sa,
                          runCol = "key",
                          quantCols = "PG.Quantity")
rm(prot)

## Add rownames

for (i in seq_along(bubis2025_prot)) {
  rownames(bubis2025_prot[[i]]) <- rowData(bubis2025_prot[[i]])$PG.ProteinGroups
  }

## Join assays
bubis2025_prot <- joinAssays(bubis2025_prot,
                  i = names(bubis2025_prot)[bubis2025_prot$fig == "fig 1"],
                  name = "proteins_fig1")
colData(bubis2025_prot[["proteins_fig1"]]) <- colData(bubis2025_prot)[bubis2025_prot$fig == "fig 1", ]

bubis2025_prot <- joinAssays(bubis2025_prot,
                       i = names(bubis2025_prot)[bubis2025_prot$fig == "fig 2"],
                       name = "proteins_fig2")
colData(bubis2025_prot[["proteins_fig2"]]) <- colData(bubis2025_prot)[bubis2025_prot$fig == "fig 2", ]

bubis2025_prot <- joinAssays(bubis2025_prot,
                       i = names(bubis2025_prot)[bubis2025_prot$fig == "fig 4"],
                       name = "proteins_fig4")
colData(bubis2025_prot[["proteins_fig4"]]) <- colData(bubis2025_prot)[bubis2025_prot$fig == "fig 4", ]

bubis2025_prot <- joinAssays(bubis2025_prot,
                       i = names(bubis2025_prot)[bubis2025_prot$fig == "fig 5"],
                       name = "proteins_fig5")
colData(bubis2025_prot[["proteins_fig5"]]) <- colData(bubis2025_prot)[bubis2025_prot$fig == "fig 5", ]

####---- Add protein data ----####

bubis2025 <- addAssay(bubis2025, bubis2025_prot[["proteins_fig1"]], "proteins_fig1")
bubis2025 <- addAssay(bubis2025, bubis2025_prot[["proteins_fig2"]], "proteins_fig2")
bubis2025 <- addAssay(bubis2025, bubis2025_prot[["proteins_fig4"]], "proteins_fig4")
bubis2025 <- addAssay(bubis2025, bubis2025_prot[["proteins_fig5"]], "proteins_fig5")

# save as rda file
save(bubis2025,
     file = file.path("../data/bubis2025.rda"),
     compress = "xz",
     compression_level = 9)


