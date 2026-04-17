library(scp)
library(tidyr)
library(dplyr)
library(ggplot2)

pep <- read.table("/mnt/disk3/leopoldguyot/PXD071075/peptide_report.tsv",
    header = TRUE,
    sep = "\t"
)
prot <- read.table("/mnt/disk3/leopoldguyot/PXD071075/protein_report.tsv",
    header = TRUE,
    sep = "\t"
)

design <- read.csv("~/dev/2025-phd-leopold-guyot/year1/scpdata_datasets/data/1.clean_meta.csv")
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
cat("PSM")
qfeatures <- readSCP(pep, runCol = "cellID", quantCols = c("Precursor.Quantity"), fnames = "Precursor.Id")
colData(qfeatures) <- as(design, "DataFrame")
psmNames <- names(qfeatures)
pepNames <- paste0("peptides_", names(qfeatures))


cat("Peptides")
qfeatures <- aggregateFeatures(qfeatures,
    i = psmNames,
    name = pepNames,
    fcol = "Modified.Sequence",
    fun = matrixStats::colMedians
)

qfeatures <- joinAssays(qfeatures, pepNames, name = "peptides")

cat("Proteins")
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

qfeatures <- addAssay(qfeatures, protse, "imported_proteins")

qfeatures <- addAssayLink(qfeatures, "peptides", "imported_proteins", varFrom = "Protein.Ids", varTo = "Protein.Ids")

saveRDS(qfeatures, "./wu2026.rds")
