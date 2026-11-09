library(QFeatures)
library(scp)
library(SummarizedExperiment)

# All data are also available on Zenodo:
# https://doi.org/0.5281/zenodo.24020222

# Update of leduc2022_plexDIA:
# https://scp.slavovlab.net/Leduc_et_al_2023

# Load PSM data
# https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?accession=MSV000093494
input_table <- read.delim(
  "report.tsv",
  header = TRUE,
  sep = "\t"
  )

# Load sample annotations
# https://drive.google.com/drive/folders/1gxBS7PlWSGzmql301hYglVgTVq1YRz-z
sample_table <- read.csv(
  "plexDIA_SingleCell_Meta_modified.csv",
  header = TRUE,
  sep = ","
  )

# Create QFeatures object
leduc2024_plexDIA <- readQFeatures(
  assayData = input_table,
  colData = sample_table,
  runCol = 'Run',
  quantCols = NULL,
  removeEmptyCols = FALSE,
  verbose = FALSE
)

leduc2024_plexDIA <- zeroIsNA(
  object = leduc2024_plexDIA,
  i = seq_along(leduc2024_plexDIA)
)

for (i in seq_along(leduc2024_plexDIA)) {
  leduc2024_plexDIA[[i]] <- logTransform(
    object = leduc2024_plexDIA[[i]],
    base = 2
  )
}

leduc2024_plexDIA <- setQFeaturesType(
  leduc2024_plexDIA,
  type = 'scp'
)

names(leduc2024_plexDIA) <- paste0(
  names(leduc2024_plexDIA),
  '_initial_import'
)

psmNames <- names(leduc2024_plexDIA)
pepNames <- paste0("peptides_", names(leduc2024_plexDIA))

# Aggregate to peptide level
leduc2024_plexDIA <- aggregateFeatures(
  leduc2024_plexDIA,
  i = psmNames,
  name = pepNames,
  fcol = "Modified.Sequence",
  fun = matrixStats::colMedians
  )

leduc2024_plexDIA <- joinAssays(
  leduc2024_plexDIA,
  pepNames,
  name = "peptides"
  )

leduc2024_plexDIA <- leduc2024_plexDIA[, , setdiff(names(leduc2024_plexDIA), pepNames)]

leduc2024_plexDIA <- addAssayLink(
  leduc2024_plexDIA,
  from = psmNames,
  to = "peptides",
  varFrom = rep("Modified.Sequence", length(psmNames)),
  varTo = "Modified.Sequence"
)

# Aggregate to protein level
leduc2024_plexDIA <- aggregateFeatures(
  leduc2024_plexDIA,
  i = "peptides",
  name = "proteins",
  fcol = "Protein.Ids",
  fun = matrixStats::colMedians
)

# Import processed protein data from paper
prot <- read.csv(
  "plexDIA_SingleCellxProtein.csv",
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
idx <- match(colnames(protse), leduc2024_plexDIA$ID)
colnames(protse) <- rownames(colData(leduc2024_plexDIA))[idx]

# Add rowData to protse
rp <- rowData(protse)
rp$Row.names <- rownames(rp)
rq <- rowData(leduc2024_plexDIA[["proteins"]])
rq$Row.names <- rq$Protein.Ids

rowData(protse) <- left_join(as.data.frame(rp), as.data.frame(rq), by = "Row.names")

# Add processed protein data from paper
leduc2024_plexDIA <- addAssay(
  leduc2024_plexDIA,
  protse,
  "imported_proteins"
  )

leduc2024_plexDIA <- addAssayLink(
  leduc2024_plexDIA,
  "peptides",
  "imported_proteins",
  varFrom = "Protein.Ids",
  varTo = "Protein.Ids"
  )

# Save QFeatures object
save(leduc2024_plexDIA,
     file = "leduc2024_plexDIA.rda",
     compress = "xz",
     compression_level = 9
)