## Ref. https://doi.org/10.1038/s41467-025-56736-7

annotations <- read.delim("./data/celltype_annotations_fresh.csv", sep = ",")
psms <- read.delim("./data/evidence.txt")
fresh <- read.csv("./data/meta_fresh.csv")

## The authors removed some of the datapoints (see the differences:)
length(unique(annotations$X)) != length(unique(fresh$ID))

## Runs are in the annotations and the identifications
sum(sort(unique(psms$Raw.file)) == sort(unique(fresh$Run)))

fresh$quantCols <- fresh$ID
annotations$quantCols <- annotations$X
metadata <- dplyr::right_join(fresh, annotations, by = "quantCols")

## Still all the runs are present
sum(sort(unique(psms$Raw.file)) == sort(unique(metadata$Run)))
metadata$runCol <- metadata$Run

## Create quantCols
quantCols <- paste0("Reporter.intensity.", 1:32)

## Create QFeatures object
leduc2025 <- scp::readSCP(psms, quantCols = quantCols, runCol = "Raw.file", colData = metadata)

## Add peptides as a summarizedExperiment and add links
peptides <- read.delim("./data/peptides.txt")
peptidesExp <- readSummarizedExperiment(peptides, quantCols = quantCols)

leduc2025 <- addAssay(leduc2025, peptidesExp, name = "Peptides")
leduc2025 <- addAssayLink(leduc2025, 1:32, "Peptides", rep("Sequence", 32), "Sequence")

## Add proteins in the same manner
proteins <- read.delim("./data/proteinGroups.txt")
proteinsExp <- readSummarizedExperiment(proteins, quantCols = quantCols)

leduc2025 <- addAssay(leduc2025, proteinsExp, name = "Proteins")
leduc2025 <- addAssayLink(leduc2025, "Peptides", "Proteins", "Proteins", "Protein.IDs")

## Save data as Rda file
save(leduc2025,
     file = "leduc2025.rda",
     compress = "xz",
     compression_level = 9)
