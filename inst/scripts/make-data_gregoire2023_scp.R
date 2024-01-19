
####---- Gregoire et al. 2023 ---####


## Grégoire, Samuel and Vanderaa, Christophe and Pyr dit Ruys,
## Sébastien and Mazzucchelli, Gabriel and Kune, Christopher and
## Vertommen, Didier and Gatto, Laurent. 2023. "Standardised workflow
## for mass spectrometry-based single-cell proteomics data processing
## and analysis using the scp package." arXiv.
## https://doi.org/10.48550/arXiv.2310.13598

## This is the script for generating the pSCoPE dataset

library(scp)
library(dplyr)
library(limma)
library(scater)

# All files were downloaded from
# https://zenodo.org/records/8417228

datadir <- "~/isilon/DataSets/SCPCBIO/MiMB_data/"

####---- Prepare sample annotations ----####

# The sample annotations are provided in separate tables for each acquisition batch:
cbio_680_annotation <- read.csv2(paste0(datadir, "sample_annotation/sampleAnnotation_cbio680.csv"))
cbio_681_annotation <- read.csv2(paste0(datadir, "sample_annotation/sampleAnnotation_cbio681.csv"))
cbio_703_annotation <- read.csv2(paste0(datadir, "sample_annotation/sampleAnnotation_cbio703.csv"))
cbio_715_annotation <- read.csv2(paste0(datadir, "sample_annotation/sampleAnnotation_cbio715.csv"))
cbio_725_annotation <- read.csv2(paste0(datadir, "sample_annotation/sampleAnnotation_cbio725.csv"))
cbio_733_annotation <- read.csv2(paste0(datadir, "sample_annotation/sampleAnnotation_cbio733.csv"))
cbio_754_annotation <- read.csv2(paste0(datadir, "sample_annotation/sampleAnnotation_cbio754bis.csv"))
giga_annotation <-  read.csv2(paste0(datadir, "sample_annotation/sampleAnnotation_GIGA.csv"))

# Bind the annotation tables into one
annotations <-
  rbind(cbio_680_annotation, cbio_681_annotation, cbio_703_annotation,
        cbio_715_annotation, cbio_725_annotation, cbio_733_annotation,
        cbio_754_annotation, giga_annotation)


####---- Prepare PSM data ----####

## Read sage id files
id_cbio <- read.delim(paste0(datadir, "sage/results.sage.cbio.tsv"))
id_giga <- read.delim(paste0(datadir, "sage/results.sage.giga.tsv"))

## Read sage quant files
quant_cbio <- read.delim(paste0(datadir, "sage/quant.cbio.tsv"))
quant_giga <- read.delim(paste0(datadir, "sage/quant.giga.tsv"))

## Bind cbio and giga data
id <- rbind(id_cbio, id_giga)
quant <- rbind(quant_cbio, quant_giga)

rm(id_cbio)
rm(id_giga)
rm(quant_cbio)
rm(quant_giga)

## Merge identification and quantification data
sage_data <- merge(quant, id,
                   by.x = c("file", "scannr"),
                   by.y = c("filename", "scannr"))

rm(quant)
rm(id)

## Add key
sage_data$.KEY <- paste(sage_data$file, sage_data$scannr, sep = ".")

## Add tmt channel name to quantitative columns
colnames(sage_data)[4:19] <- paste0(c(126, rep(127:133, each = 2), 134),
                                    rep(c("C", "N"), 8))

## Add run variable
sage_data$run <- sub(".*_1_1_", "", sage_data$file)
sage_data$run <- sub("^1", "GIGA_1", sage_data$run)
sage_data$run <- sub("\\.mzML", "", sage_data$run)

## Add experiment batches
sage_data$batch <- sub("_.*", "", sage_data$run)

gregoire2023_scp <- readSCP(featureData = sage_data,
                            colData = annotations,
                            channelCol = "channel",
                            batchCol = "run",
                            removeEmptyCols = TRUE,
                            sep = "_")

####---- Retrieve processed data ----####

# Replace 0s by NAs
gregoire2023_scp <- zeroIsNA(gregoire2023_scp, i = 1:length(gregoire2023_scp))

# Keep targets, ranks 1 and fdr < 0.01
gregoire2023_scp <- filterFeatures(gregoire2023_scp,
                      ~ rank == 1 &
                        peptide_fdr < 0.01 &
                        label == 1)

# Highlight chimeric spectra
for (i in seq_along(gregoire2023_scp)) {
  # Extract rowData for each set
  rd <- rowData(gregoire2023_scp[[names(gregoire2023_scp)[i]]])
  # Create unique spectrum identifier .KEY
  rd$.KEY <- paste(rd$file, rd$scannr)
  # Create "chimeric" column, FALSE by default
  rd$chimeric <- FALSE
  # Change "chimeric" to TRUE for duplicated keys
  rd$chimeric[rd$.KEY %in% rd$.KEY[duplicated(rd$.KEY)]] <- TRUE
  # Store updated rowData
  rowData(gregoire2023_scp[[names(gregoire2023_scp)[i]]]) <- rd
}

# Remove chimeric spectra
gregoire2023_scp <- filterFeatures(gregoire2023_scp,
                      ~ !chimeric)

# Compute sample to carrier ratio
gregoire2023_scp <- computeSCR(gregoire2023_scp,
                  i = 1:length(gregoire2023_scp),
                  colvar = "cell_type",
                  carrierPattern = "carrier",
                  samplePattern = "THP1|THP1_dif|U937|U937_dif|mix",
                  rowDataName = "MeanSCR")

# Remove SCR > 1
gregoire2023_scp <- filterFeatures(gregoire2023_scp,
                      ~ !is.na(MeanSCR) &
                        MeanSCR < 1)

# Remove carrier and empty channel
gregoire2023_scp <- subsetByColData(gregoire2023_scp, !gregoire2023_scp$cell_type %in% c("carrier", "empty"))

# Compute median intensity
for (i in names(gregoire2023_scp)) {
  # Extract log assay
  logAssay <- log(assay(gregoire2023_scp[[i]]))
  # Compute median RI by cell
  meds <- colMedians(logAssay, na.rm = TRUE, useNames = TRUE)
  # Store median RI in colData.
  colData(gregoire2023_scp)[names(meds), "log_medianRI"] <- meds
}

# Compute median coefficient of variation
gregoire2023_scp <- medianCVperCell(gregoire2023_scp,
                       i = 1:length(gregoire2023_scp),
                       groupBy = "proteins",
                       nobs = 3,
                       norm = "div.median",
                       colDataName = "medianCV")

# Count number of peptide per cell
gregoire2023_scp <- countUniqueFeatures(gregoire2023_scp,
                           i = 1:length(gregoire2023_scp),
                           groupBy = "peptide",
                           colDataName = "count")

# filter poor quality cells
filter_samples <-
  (gregoire2023_scp$batch == "CBIO680" & gregoire2023_scp$log_medianRI > 7.77 &
     gregoire2023_scp$count > 1250 & gregoire2023_scp$medianCV < 0.615) |
  (gregoire2023_scp$batch == "CBIO681" & gregoire2023_scp$log_medianRI > 8.5 &
     gregoire2023_scp$count > 1250 & gregoire2023_scp$medianCV < 0.79) |
  (gregoire2023_scp$batch == "CBIO703" & gregoire2023_scp$log_medianRI > 7.69 &
     gregoire2023_scp$count > 1250 & gregoire2023_scp$medianCV < 0.68) |
  (gregoire2023_scp$batch == "CBIO715" & gregoire2023_scp$log_medianRI > 7.69 &
     gregoire2023_scp$count > 1250 & gregoire2023_scp$medianCV < 0.62) |
  (gregoire2023_scp$batch == "CBIO725" & gregoire2023_scp$log_medianRI > 8.08 &
     gregoire2023_scp$count > 1250 & gregoire2023_scp$medianCV < 0.73) |
  (gregoire2023_scp$batch == "CBIO754" & gregoire2023_scp$log_medianRI > 7.39 &
     gregoire2023_scp$count > 1250 & gregoire2023_scp$medianCV < 0.67) |
  (gregoire2023_scp$batch == "GIGA" &
     gregoire2023_scp$count > 1250 & gregoire2023_scp$medianCV < 0.455)

gregoire2023_scp <- subsetByColData(gregoire2023_scp, filter_samples) |>
  dropEmptyAssays()

# Rmove blanks
gregoire2023_scp <- subsetByColData(gregoire2023_scp, gregoire2023_scp$cell_type != "blank")

# Aggregate PSMs into peptides
gregoire2023_scp <- aggregateFeatures(gregoire2023_scp,
                         i = 1:length(gregoire2023_scp),
                         fcol = "peptide",
                         name = paste0("peptide_", names(gregoire2023_scp)),
                         fun = colMedians, na.rm = TRUE)

# Join PSM assays per acquisition batches
batches <- c("CBIO680", "CBIO681", "CBIO703",
             "CBIO715", "CBIO725", "CBIO754",
             "GIGA")

for (batch in batches) {
  gregoire2023_scp <- joinAssays(gregoire2023_scp,
                    i = grep(paste0("peptide_", batch), names(gregoire2023_scp)),
                    name = paste0("peptides_", batch))
}

# Remove highly missing peptides
gregoire2023_scp <- filterNA(gregoire2023_scp,
                i = grep("peptides", names(gregoire2023_scp)),
                pNA = 0.98)

pep_assay_names <- names(gregoire2023_scp)[grep("peptides_", names(gregoire2023_scp))]

# Normalise cells
for (i in seq_along(pep_assay_names)) {
  gregoire2023_scp <- sweep(gregoire2023_scp,
               i = pep_assay_names[i],
               MARGIN = 2,
               FUN = "/",
               STATS = colMedians(assay(gregoire2023_scp[[pep_assay_names[i]]]), na.rm = TRUE),
               name = paste0(pep_assay_names[i], "_norm"))
}

# Log transform normalised peptide assays
pep_assay_names <- names(gregoire2023_scp)[grep("peptides_.*_norm", names(gregoire2023_scp))]

gregoire2023_scp <- logTransform(gregoire2023_scp,
                    base = 2,
                    i = pep_assay_names,
                    name = paste0(pep_assay_names, "_log"))

pep_assay_names <- names(gregoire2023_scp)[grep("peptides_.*_norm_log", names(gregoire2023_scp))]

gregoire2023_scp <- aggregateFeatures(gregoire2023_scp,
                         i = pep_assay_names,
                         fcol = "proteins",
                         fun = colMedians, na.rm = TRUE,
                         name = sub("peptides", "proteins", pep_assay_names))

# Remove batch effects
for (i in grep("norm_log|imptd", names(gregoire2023_scp))) {
  ## Extract set
  sce <- getWithColData(gregoire2023_scp, names(gregoire2023_scp)[i])
  ## Batch correct assay
  assay(sce) <-
    removeBatchEffect(assay(sce), group = sce$cell_type,
                      batch = sce$run, batch2 = sce$channel)
  ## Name and add batch-corrected assay
  gregoire2023_scp <- addAssay(gregoire2023_scp,
                  y = sce,
                  name = sub("_norm_log|mptd", "_batchC", names(gregoire2023_scp)[i]))
  ## Add link between batch corrected and original assay
  gregoire2023_scp <- addAssayLinkOneToOne(gregoire2023_scp,
                              from = names(gregoire2023_scp)[i],
                              to = sub("_norm_log|mptd", "_batchC", names(gregoire2023_scp)[i]))
}

# Compute NIPALS dimensionality reduction
for (i in grep("batchC", names(gregoire2023_scp))) {
  nipals_res <-
    ## Extract assay
    assay(gregoire2023_scp[[i]]) |>
    as.data.frame() |>
    ## Encode missing values
    mutate_all(~ifelse(is.nan(.), NA, .)) |>
    ## Transpose
    t() |>
    ## PCA
    pcaMethods::pca(method="nipals", nPcs = 2)

  reducedDim(gregoire2023_scp[[i]], "NIPALS") <- pcaMethods::scores(nipals_res)
}

# Save data as Rda file
save(gregoire2023_scp,
     file = "~/2022-phd-samuel-gregoire/gregoire2023_scp.Rda",
     compress = "xz",
     compression_level = 9)
