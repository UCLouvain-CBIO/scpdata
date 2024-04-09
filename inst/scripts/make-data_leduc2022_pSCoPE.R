
####---- Leduc et al. 2022 ---####


## Leduc, Andrew, R. Gray Huffman, and Nikolai Slavov. 2021. “Droplet
## Sample Preparation for Single-Cell Proteomics Applied to the Cell
## Cycle.” bioRxiv. https://doi.org/10.1101/2021.04.24.441211.

## This is the script for generating the pSCoPE dataset

library("SingleCellExperiment")
library("scp")
library("tidyverse")


# All files were downloaded from
# https://drive.google.com/drive/folders/117ZUG5aFIJt0vrqIxpKXQJorNtekO-BV

datadir <- "~/Documents/.localData/SCP/leduc2022/pSCoPE/"

####---- Prepare sample annotations ----####

# Get LC and MS annotations with sample type annotations
design <- read.csv(paste0(datadir, "CellenONE_data/annotation_18plex_well.csv"))
batch <- read.csv(paste0(datadir, "batch.csv"))

# Clean the sample metadata so that it meets the requirements for
# `scp::readSCP`. We first need to transform the design (set x
# reporter ion) to a long table so that one line is one sample.
design <- pivot_longer(design, -c(Set, well), names_to = "Channel",
                       values_to = "SampleAnnotation")
design$SampleType <- recode(design$SampleAnnotation,
                            neg = "Negative",
                            u = "Monocyte",
                            m = "Melanoma",
                            unused = "Unused",
                            reference = "Reference",
                            carrier = "Carrier")

# We then make some slight corrections to the batch data
colnames(batch)[1] <- "Set" ## consistent naming with design
batch$digest <- as.character(batch$digest)

# We can now combine the two tables in a single annotation table
sampleAnnotation <- inner_join(design, batch, by = "Set")

# Add melanoma subpopulation information
index <- read.csv(paste0(datadir, "misc/sample_index.csv"), row.names = 1)
meta <- read.csv(paste0(datadir, "misc/meta.csv"), row.names = 1)
meta <- inner_join(index, meta, by = "id")

idAnnot <- paste0(sampleAnnotation$Set, sub("RI", "", sampleAnnotation$Channel))
idMeta <- paste0(sub("^X", "", meta$rawfile), meta$channel.x)

sampleAnnotation$MelanomaSubCluster <- meta$sub[match(idAnnot, idMeta)]
sampleAnnotation$MelanomaSubCluster <- recode(sampleAnnotation$MelanomaSubCluster, C1 = "A", C2 = "B")

# Add cellenONE information = sample prep annotations
# 2 utility functions:
parseCellenoneFieldData <- function(file, ...) {
    out <- read.table(
        file, sep = "\t", fill = TRUE, header = FALSE,
        col.names = c("position", "well","volume"),  quote = "", ...
    )
    containsFieldData <- grepl("\\[\\d", out$position)
    fieldData <- out$position[containsFieldData]
    out$Field[containsFieldData] <- sub("\\[(\\d*),.*$", "\\1", fieldData)
    out$Field <- as.numeric(out$Field) + 1
    out <- fill(out, Field, .direction = "down")
    out <- out[out$well != "", ]
    out$YPos <- as.numeric(sub("^(.*)/.*$", "\\1", out$position))
    out$XPos <- as.numeric(sub("^.*/(.*)$", "\\1", out$position))
    out
}
matchLayoutToPickup <- function(layout, pickup,
                                coords = c("XPos","YPos")) {
    require("yaImpute")
    out <- lapply(split(layout, layout$Field), function(x) {
        field <- unique(x$Field)
        ref <- unique(pickup[pickup$Field == field, coords])
        target <- x[, coords]
        ann <- ann(
            ref = as.matrix(ref), target = as.matrix(target), k = 1
        )
        pickupIndex <- ann$knnIndexDist[, 1]
        x$XPosPickup <- ref[pickupIndex,]$XPos
        x$YPosPickup <- ref[pickupIndex,]$YPos
        x$PickupIndex <- as.character(pickupIndex)
        x
    })
    do.call(rbind, out)
}

# Parse cell isolation files
macrophage <- read.table(
    file = paste0(datadir, "CellenONE_data/M_isolated.csv"),
    sep = ",", header = TRUE
)
macrophage$Condition <- "macrophage"
monocyte <- read.table(
    file = paste0(datadir, "CellenONE_data/U_isolated.csv"),
    sep = ",", header = TRUE
)
monocyte$Condition <- "monocyte"
cellenoneData <- rbind(macrophage, monocyte)
isOdd <- 1:nrow(cellenoneData) %% 2 == 1
cellenoneData <- cbind(
    cellenoneData[!isOdd, 1:7],
    cellenoneData[isOdd, 8:ncol(cellenoneData)]
)
cellenoneData$IsolationTimeStamp <- paste(cellenoneData$Date, cellenoneData$Time)
cellenoneData$IsolationTimeStamp <- sub(
    "/21[ ]", "/2021 ", cellenoneData$IsolationTimeStamp
)
cellenoneData$IsolationTimeStamp <- strptime(
    cellenoneData$IsolationTimeStamp, "%m/%d/%Y %I:%M:%S %p"
)

# Parse labelling information
label <- parseCellenoneFieldData(
    paste0(datadir, "CellenONE_data/Labels.txt"), skip = 22
)
label$Label <- as.numeric(sub("^.*[A-Z](.*),", "\\1", label$well)) + 4
label$Label <- paste0("RI", label$Label)
label <- dplyr::rename(label, WellLabel = well)
label$WellLabel <- gsub("^1|,", "", label$WellLabel)
cellenoneData <- left_join(
    cellenoneData, label, by = join_by(XPos, YPos, Field)
)

# Parse sample pick-up for LC-MS information
pickup <- parseCellenoneFieldData(
    paste0(datadir, "CellenONE_data/SamplePickup.txt"), skip = 27
)
pickup$Target <- 5 - ceiling(pickup$Field / 4)
pickup$Field <- rep(c(1,2,3,4), each = 9)
pickup$well <- sub("^.*([A-Z].*),", "\\1", pickup$well)
pickup <- dplyr::rename(pickup, WellPooled = well)
cellenoneData <- matchLayoutToPickup(cellenoneData, pickup)
cellenoneData <- left_join(
    cellenoneData, pickup,
    join_by(XPosPickup == XPos, YPosPickup == YPos, Field, Target)
)
cellenoneData <- dplyr::select(
    cellenoneData,
    IsolationTimeStamp, Diameter, Elongation ,Target, Field, XPos, YPos,
    XPosPickup, YPosPickup, WellPooled, Label
)
cellenoneData <- dplyr::rename(
    cellenoneData, GlassSlide = Target, XPosDrop = XPos, YPosDrop = YPos,
    Channel = Label,
)

# Combine all sample annotations
sampleAnnotation <- left_join(
    sampleAnnotation, cellenoneData,
    by = join_by(well == WellPooled, Channel)
)
sampleAnnotation <- dplyr::rename(sampleAnnotation, WellPooled = well)

sampleAnnotation |>
    filter(!is.na(Field) & !is.na(GlassSlide)) |>
    ggplot() +
    ggforce::geom_circle(aes(x0 = XPosPickup,
                    y0 = YPosPickup,
                    fill = lcbatch,
                    r = 10),
                color = "transparent") +
    geom_point(aes(x = XPosDrop,
                   y = YPosDrop,
                   colour = SampleType
                   )) +
    geom_text(aes(x = XPosPickup,
                  y = YPosPickup,
                  label = WellPooled),
              size = 3, colour = "white") +
    facet_grid(GlassSlide ~ Field, labeller = label_both) +
    theme_minimal() +
    scale_fill_manual(values = c("palegreen", "lightgoldenrod"))

####---- Prepare PSM data ----####

ev <- read.delim(paste0(datadir, "ev_updated.txt"))
colnames(ev) <- gsub("^Reporter.intensity.(\\d*)$", "RI\\1", colnames(ev))
colnames(ev)[colnames(ev) == "Raw.file"] <- "Set"
ev$modseq <- paste0(ev$Modified.sequence, ev$Charge)
## This removes DIA runs
ev <- ev[ev$Set %in% sampleAnnotation$Set, ]

## Create the QFeatures object
leduc2022_pSCoPE <- readSCP(ev, sampleAnnotation,
                            channelCol = "Channel",
                            batchCol = "Set")

## Clean protein names
rdList <- lapply(rowData(leduc2022_pSCoPE), function(rd) {
    rd$Leading.razor.protein.id <-
        gsub("^.*\\|(.*)\\|.*", "\\1", rd$Leading.razor.protein)
    rd$Leading.razor.protein.symbol <-
        gsub("^.*\\|.*\\|(.*)_.*", "\\1", rd$Leading.razor.protein)
    rd
})
rowData(leduc2022_pSCoPE) <- rdList

####---- Retrieve processed data ----####

## Retrieve the data processed by Leduc et al.
sampleInd <- read.csv(paste0(datadir, "misc/sample_index.csv"), row.names = 2)
files <- c("t0.csv", "t3.csv", "t4b.csv", "t6.csv")
processedData <- lapply(files, function(f) {
    ## Read data
    dat <- read.csv(paste0(datadir, "processed_data/", f), row.names = 1)
    dat <- as.matrix(dat)
    ## Convert column names
    fileID <- sub("X", "", sampleInd[colnames(dat), "rawfile"])
    channel <- sampleInd[colnames(dat), "channel"]
    colnames(dat) <- paste0(fileID, "RI", channel)
    ## Convert to a SCE
    dat <- SingleCellExperiment(dat)
    ## Add colData
    colData(dat) <- colData(leduc2022_pSCoPE)[colnames(dat), ]
    dat
})
names(processedData) <- c("peptides", "peptides_log", "proteins_norm2", "proteins_processed")

## Generate the peptide to protein table
pep2prot <- rbindRowData(leduc2022_pSCoPE, names(leduc2022_pSCoPE)) %>%
    data.frame %>%
    group_by(modseq) %>%
    summarise(Leading.razor.protein = paste(unique(Leading.razor.protein),
                                            collapse = ";"),
              Leading.razor.protein.id = paste(unique(Leading.razor.protein.id),
                                               collapse = ";"),
              Leading.razor.protein.symbol = paste(unique(Leading.razor.protein.symbol),
                                                   collapse = ";")) %>%
    data.frame
rownames(pep2prot) <- pep2prot$modseq

## Add `peptides` data
rowData(processedData$peptides) <- pep2prot[rownames(processedData$peptides), ]
leduc2022_pSCoPE <- addAssay(leduc2022_pSCoPE, processedData$peptides, name = "peptides")
leduc2022_pSCoPE <- addAssayLink(leduc2022_pSCoPE, from = 1:134, to = "peptides",
                      varFrom = rep("modseq", 134), varTo = "modseq")

## Add `peptides_log` data
rowData(processedData$peptides_log) <- pep2prot[rownames(processedData$peptides_log), ]
leduc2022_pSCoPE <- addAssay(leduc2022_pSCoPE, processedData$peptides_log, name = "peptides_log")
leduc2022_pSCoPE <- addAssayLink(leduc2022_pSCoPE, from = "peptides", to = "peptides_log",
                      varFrom = "modseq", varTo = "modseq")

## Add `proteins_norm2` data
prots <- select(pep2prot, Leading.razor.protein, Leading.razor.protein.id, Leading.razor.protein.symbol)
prots <- prots[!duplicated(prots$Leading.razor.protein.id), ]
rownames(prots) <- prots$Leading.razor.protein.id
rowData(processedData$proteins_norm2) <- prots[rownames(processedData$proteins_norm2), ]
leduc2022_pSCoPE <- addAssay(leduc2022_pSCoPE, processedData$proteins_norm2, name = "proteins_norm2")
leduc2022_pSCoPE <- addAssayLink(leduc2022_pSCoPE, from = "peptides_log", to = "proteins_norm2",
                      varFrom = "Leading.razor.protein.id",
                      varTo = "Leading.razor.protein.id")

## Add `proteins_processed` data
rowData(processedData$proteins_processed) <- prots[rownames(processedData$proteins_processed), ]
leduc2022_pSCoPE <- addAssay(leduc2022_pSCoPE, processedData$proteins_processed, name = "proteins_processed")
leduc2022_pSCoPE <- addAssayLinkOneToOne(leduc2022_pSCoPE, from = "proteins_norm2", to = "proteins_processed")

## Create more informative TMT labels
tmtlabs <- c(
    "TMT126", paste0(
        "TMT", rep(127:134, each = 2), rep(c("N", "C"), 8)
    ), "TMT135N"
)
leduc2022_pSCoPE$Channel <- factor(
    tmtlabs[as.numeric(sub("RI", "", leduc2022_pSCoPE$Channel))],
    levels = tmtlabs
)

# Save data as Rda file
save(leduc2022_pSCoPE,
     file = "~/Documents/.localData/scpdata/leduc2022_pSCoPE.Rda",
     compress = "xz",
     compression_level = 9)
