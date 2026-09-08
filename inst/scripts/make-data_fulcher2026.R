
####---- Fulcher et al. 2026 ---####

library("SingleCellExperiment")
library("scp")
library("dplyr")
library("tidyr")

## TODO replace with Zenodo link and use BiocFileCache
dataDir <- "~/Documents/research/.localData/SCP/fulcher2026/"

####---- Quantitative data ----####

peptides <- paste0(dataDir, "msv000100684/FragPipe/tmt-report/abundance_peptide_None.tsv") |>
    read.delim(, check.names = FALSE)
proteins <- paste0(dataDir, "msv000100684/FragPipe/tmt-report/abundance_protein_None.tsv") |>
    read.delim(, check.names = FALSE)

####---- Sample annotations ----####

## Get sample names
annots <- data.frame(
    SampleName = grep("PBMCs", colnames(proteins), value = TRUE)
) |>
    ## Extract annotations from sample names
    mutate(
        PoolId = sub("_\\d.+?$", "", SampleName),
        Well = sub(".*_([A-Z]\\d+)_.*", "\\1", SampleName),
        Channel = sub(".*_(.*)$", "\\1", SampleName),
        Chip = case_when(
            grepl("Target3", SampleName) ~ "Target3",
            grepl("Target4_A[12]_", SampleName) ~ "Target3",
            .default = "Target4"
        ),
        ChannelType = case_when(
            Channel %in% c("127C", "134ND") ~ "Isotope",
            Channel %in% c("135CD", "134C", "135N") ~ "Blank",
            grepl("^Ref", SampleName) ~ "Bridge",
            .default = "Single-Cell"
        )
    ) |>
    ## Add LC batch information
    left_join(
        paste0(dataDir, "github/Single_PBMCs_Fig4/Metadata/LC_Column_Meta.csv") |>
            read.csv() |>
            rename(PoolId = File_Name),
        by = "PoolId"
    )
## Add the cell filter derived by the authors
kept <- paste0(dataDir, "github/Single_PBMCs_Fig4/Metadata/First_Pass_PBMCs_KEEP.csv") |>
    read.csv() |>
    pull(SampleID)
annots$RemovedByAuthors <- !annots$SampleName %in% kept

## There is no data to map the cellenone data to the ms data...
isolation <- read.delim(paste0(dataDir, "zenodo_18421469/20241218_130552.609__isolated.xls"))
nevents <- nrow(isolation)
isolation <- cbind(
    isolation[seq(1, nevents, 2), c(1, 8:16)],
    isolation[seq(2, nevents, 2), 2:7]
)

####---- Format to a QFeatures ----####

## Read peptide data
annots$quantCols <- annots$SampleName
fulcher2026 <- readSCP(
    peptides, colData = annots, name = "peptides", fnames = "Peptide"
)

## Add protein data
se <- readSummarizedExperiment(
    proteins, quantCols = annots$quantCols, fnames = "Protein"
)
fulcher2026 <- addAssay(fulcher2026, se, "proteins")
fulcher2026 <- addAssayLink(
    fulcher2026,
    from = "peptides",
    to = "proteins",
    varFrom = "Protein",
    varTo = "Protein"
)

####---- Test data with scplainer ----####

library("ggplot2")
library("patchwork")

## Minimal processing
se <- getWithColData(fulcher2026,  "peptides")
se <- zeroIsNA(se)
se <- logTransform(se)
se$NPeptides <- colSums(!is.na(assay(se)))
se$NProteins <- colSums(!is.na(assay(se)))
se$MedianIntensity <- colMedians(assay(se), na.rm = TRUE)
colData(se) |>
    data.frame() |>
    ggplot() +
    aes(
        x = NPeptides,
        y = MedianIntensity,
        colour = ChannelType,
        shape = RemovedByAuthors
    ) +
    geom_point()
se <- se[, !se$RemovedByAuthors & ChannelType == "Single-Cell"]

## scplainer
se <- scpModelWorkflow(
    se,
    formula = ~ 1 + ## intercept
        MedianIntensity + ## normalization
        ## batch effects
        Channel + PoolId
)
scpModelFilterPlot(se)
scpModelFilterThreshold(se) <- 5
(caRes <- scpComponentAnalysis(
    se, ncomp = 2, method = "APCA", effect = character(0)
))
caResCells <- caRes$bySample
se$cell <- colnames(se)
caResCells <- scpAnnotateResults(caResCells, colData(se), by = "cell")
## Plot results
scpComponentPlot(
    caResCells,
    pointParams = list(aes(colour = NPeptides, shape = as.factor(LC)))
) |>
    wrap_plots() +
    plot_layout(guides = "collect") &
    scale_colour_continuous(type = "viridis")
## Add modelled data back to QFeatures object
## Keep only set-specific cell annotations
keep <- c("NPeptides", "NProteins", "MedianIntensity")
colData(se) <- colData(se)[, keep]
fulcher2026 <- addAssay(fulcher2026, se, "proteins_modelled")
fulcher2026 <- addAssayLink(
    fulcher2026,
    from = "proteins",
    to = "proteins_modelled",
    varFrom = "Protein",
    varTo = "Protein"
)

####---- Save data ----####

# Save data as Rda file
save(
     fulcher2026,
     file = "~/Documents/research/.localData/scpdata/fulcher2026.rda",
     compress = "xz",
     compression_level = 9
)


