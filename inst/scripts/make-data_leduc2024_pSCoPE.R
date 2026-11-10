library(QFeatures)
library(scp)
library(SummarizedExperiment)
library(fuzzyjoin)

# All data are also available on Zenodo:
# https://doi.org/10.5281/zenodo.24020112

# Update of leduc2022_pSCoPE:
# https://scp.slavovlab.net/Leduc_et_al_2023

# Load PSM data
# https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?accession=MSV000094207
input_table <- read.delim(
  "evidence.txt",
  header = TRUE,
  sep = "\t"
  )

# Load sample annotations
# https://drive.google.com/drive/folders/10LYxtc7qNdO06z8FKtf0L9UB4IAta68h
sample_table <- read.csv(
  "TMT_SingleCell_Meta_modified.csv",
  header = TRUE,
  sep = ","
  )

# Create QFeatures object
leduc2024_pSCoPE <- readQFeatures(
  assayData = input_table,
  colData = sample_table,
  runCol = 'Raw.file',
  quantCols = NULL,
  removeEmptyCols = FALSE,
  verbose = FALSE
)

leduc2024_pSCoPE <- zeroIsNA(
  object = leduc2024_pSCoPE,
  i = seq_along(leduc2024_pSCoPE)
)

for (i in seq_along(leduc2024_pSCoPE)) {
  leduc2024_pSCoPE[[i]] <- logTransform(
    object = leduc2024_pSCoPE[[i]],
    base = 2
  )
}

leduc2024_pSCoPE <- setQFeaturesType(
  leduc2024_pSCoPE,
  type = 'scp'
)

names(leduc2024_pSCoPE) <- paste0(
  names(leduc2024_pSCoPE),
  '_initial_import'
)

psmNames <- names(leduc2024_pSCoPE)
pepNames <- paste0("peptides_", names(leduc2024_pSCoPE))

# Aggregate to peptide level
leduc2024_pSCoPE <- aggregateFeatures(
  leduc2024_pSCoPE,
  i = psmNames,
  name = pepNames,
  fcol = "Modified.sequence",
  fun = matrixStats::colMedians
)

leduc2024_pSCoPE <- joinAssays(
  leduc2024_pSCoPE,
  pepNames, name = "peptides"
  )

leduc2024_pSCoPE <- leduc2024_pSCoPE[, , setdiff(names(leduc2024_pSCoPE), pepNames)]

leduc2024_pSCoPE <- addAssayLink(
  leduc2024_pSCoPE,
  from = psmNames,
  to = "peptides",
  varFrom = rep("Modified.sequence", length(psmNames)),
  varTo = "Modified.sequence"
)

# Aggregate to protein level
leduc2024_pSCoPE <- aggregateFeatures(
  leduc2024_pSCoPE,
  i = "peptides",
  name = "proteins",
  fcol = "Leading.razor.protein",
  fun = matrixStats::colMedians
)

# Import processed protein data from paper
prot <- read.csv(
  "TMT_ProteinxSingleCell.csv",
  header = TRUE,
  sep = ",",
  row.names = 1
  )

# Create SummarizedExperiment object
counts <- as.matrix(prot)
protse <- SummarizedExperiment(
  assays = list(counts = counts)
)

# Match column names in quantification and sample annotation tables

# Match protse colnames with QF sample names
idx <- match(colnames(protse), leduc2024_pSCoPE$ID)
colnames(protse) <- rownames(colData(leduc2024_pSCoPE))[idx]

# Add rowData to protse
rp <- rowData(protse)
rp$Row.names <- rownames(rp)
rq <- rowData(leduc2024_pSCoPE[["proteins"]])
rq$Row.names <- rownames(rq)

joined <- fuzzyjoin::regex_left_join(as.data.frame(rp), as.data.frame(rq), by = "Row.names")
rowData(protse) <- joined[!duplicated(joined$Row.names.x),]

# Add processed protein data from paper
leduc2024_pSCoPE <- addAssay(
  leduc2024_pSCoPE,
  protse,
  "imported_proteins"
  )

leduc2024_pSCoPE <- addAssayLink(
  leduc2024_pSCoPE,
  "peptides", 
  "imported_proteins",
  varFrom = "Leading.razor.protein",
  varTo = "Leading.razor.protein"
  )

# Save QFeatures object
save(leduc2024_pSCoPE,
     file = "leduc2024_pSCoPE.rda",
     compress = "xz",
     compression_level = 9
)
