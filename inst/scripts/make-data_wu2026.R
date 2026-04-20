library(scp)
library(tidyr)
library(dplyr)

pep <- read.table("/mnt/disk3/leopoldguyot/PXD071075/peptide_report.tsv",
    header = TRUE,
    sep = "\t"
)
prot <- read.table("/mnt/disk3/leopoldguyot/PXD071075/protein_report.tsv",
    header = TRUE,
    sep = "\t"
)

design <- read.csv("~/dev/2025-phd-leopold-guyot/year1/scpdata_datasets/data/1.clean_meta.csv")

# Prepare peptide table
rownames(design) <- design$ID
pep <- pep %>%
    separate(Run,
        into = c("Platform", "Sample", "Number", "Unknown"),
        sep = "_",
        remove = FALSE
    ) %>%
    mutate(Precursor.Quantity = log2(Precursor.Quantity))

pep$Sample <- dplyr::recode(
    pep$Sample,
    "WLBSC"  = "S1",
    "WLBSC2" = "S2",
    "WLBSC3" = "S3"
)

pep$cellID <- paste(pep$Sample, pep$Number, sep = "_")

# Loading PSM
wu2026 <- readSCP(pep, runCol = "cellID", quantCols = c("Precursor.Quantity"), fnames = "Precursor.Id")
# add colData
colData(wu2026) <- as(design, "DataFrame")

psmNames <- names(wu2026)
pepNames <- paste0("peptides_", names(wu2026))


# Aggregate to peptide level
wu2026 <- aggregateFeatures(wu2026,
    i = psmNames,
    name = pepNames,
    fcol = "Modified.Sequence",
    fun = matrixStats::colMedians
)

wu2026 <- joinAssays(wu2026, pepNames, name = "peptides")
wu2026 <- wu2026[, , setdiff(names(wu2026), pepNames)]
wu2026 <- addAssayLink(wu2026,
    from = psmNames,
    to = "peptides",
    varFrom = rep("Modified.Sequence", length(psmNames)),
    varTo = "Modified.Sequence"
)

# Aggregate to protein level
wu2026 <- aggregateFeatures(wu2026,
    i = "peptides",
    name = "proteins",
    fcol = "Protein.Ids",
    fun = matrixStats::colMedians
)

# Add paper's protein data
protse <- readSummarizedExperiment(prot,
    quantCols = grep("Evosep", colnames(prot)),
    fnames = "Protein.Ids"
)

old_names <- colnames(protse)

mapping <- c(
    "WLBSC" = "S1",
    "WLBSC2" = "S2",
    "WLBSC3" = "S3"
)

rename_col <- function(name) {
    wbsc_type <- grep(paste(names(mapping), collapse = "|"), name, value = TRUE)
    wbsc_type <- sub(".*(WLBSC[23]?).*", "\\1", name) # extract WLBSC/WLBSC2/WLBSC3
    prefix <- mapping[wbsc_type]

    # extract first number after WLBSC
    num <- sub(".*WLBSC[23]?_([0-9]+)_.*", "\\1", name)

    paste0(prefix, "_", num)
}

colnames(protse) <- sapply(old_names, rename_col)

protse <- protse[!duplicated(rownames(protse))]

wu2026 <- addAssay(wu2026, protse, "imported_proteins")

wu2026 <- addAssayLink(wu2026, "peptides", "imported_proteins", varFrom = "Protein.Ids", varTo = "Protein.Ids")

save(wu2026,
    file = "/mnt/disk3/leopoldguyot/wu2026.rda",
    compress = "xz",
    compression_level = 9
)
