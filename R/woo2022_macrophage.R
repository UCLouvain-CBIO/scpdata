##' Woo et al. 2022 (Cell Syst.): LPS-treated macrophages
##'
##' Single-cell data from macrophages subjected to 3 LPS
##' treatments. The data were
##' acquired using the TIFF (transfer identification based on FAIMS
##' filtering) acquisition method. The data contain 155 single cells:
##' 54 control cells (no treatment), 52 cells treated with LPS during
##' 24h and 49 cells treated with LPS during 49h.
##'
##' @format A [QFeatures] object with 5 assays, each assay being a
##' [SingleCellExperiment] object:
##'
##' - `peptides_[intensity or LFQ]`: 2 assays containing peptide
##'   quantities or normalized quantities using the maxLFQ method
##'   as computed by MaxQuant.
##' - `proteins_[intensity or iBAQ or LFQ]`: 3 assays containing
##'   protein quantities or normalized proteins using the iBAQ or
##'   maxLFQ methods as computed by MaxQuant.
##'
##' Sample annotation is stored in `colData(woo_macrophage())`.
##'
##' @section Acquisition protocol:
##'
##' The data were acquired using the following setup. More information
##' can be found in the source article (see `References`).
##'
##' - **Sample isolation**: cultured RAW 264.7 cells treated or not
##'   with 100 ng/ul LPS. The cells were sorted using the Influx II
##'   cell sorter and deposited on a nanoPOTS chip.
##' - **Sample preparation**: cells are lysed using using a DDM+DTT
##'   lysis and reduction buffer. The proteins are alkylated with IAA
##'   and digested with LysC and trypsin. Samples are then acidified
##'   with FA, vacuum dried and stored in freezer until data
##'   acquisition.
##' - **Liquid chromatography**: peptides are loaded using an in-house
##'   autosampler (Williams et al. 2020). The samples are concentrated
##'   through a SPE column (4cm x 100µm i.d. packed with 5µm C18) with
##'   microflow LC pump. The peptides are then eluted from an LC
##'   column (25cm x 50 µm i.d. packed with 1.7µm C18) from a 60 min
##'   gradient (100nL/min).
##' - **Mass spectrometry**: MS/MS was performed on an Orbitrap Fusion
##'   Lumos Tribrid MS with FAIMSpro coupled to a 2.4 kV ESI. FAIMS
##'   setup: 4-CV method (-45, -55, -65, -75 V). MS1 setup: resolution
##'   = 120.000, range = 350-1500 m/z,AGC target of 1E6, accumulation
##'   of 254ms. MS2 setup: 30% HCD, resolution AGC 2E4, accumulation
##'   of 254ms.
##' - **Raw data processing**: preprocessing using Maxquant v1.6.2.10
##'   that use Andromeda search engine (with UniProtKB 2016-21-29).
##'   MBR was enabled.
##'
##' @section Data collection:
##'
##' All data were collected from the MASSIVE repository (accession ID:
##' MSV000085937).
##'
##' The peptide and protein data were extracted from the
##' `peptides_RAW_LPS_scProteomics.txt` or
##' `proteinGroups_RAW_LPS_scProteomics.txt` files, respectively, in
##' the `RAW_LPS_SingleCellProteomics` folders.
##'
##' The tables were split so that intensities, maxLFQ, and iBAQ
##' data are contained in separate tables. Tables are then
##' converted to [SingleCellExperiment] objects. Sample annotations
##' were inferred from the sample names. All data is combined in
##' a [QFeatures] object. [AssayLinks] were stored between peptide
##' assays and their corresponding proteins assays based on the
##' leading razor protein (hence only unique peptides are linked to
##' proteins).
##'
##' The script to reproduce the `QFeatures` object is available at
##' `system.file("scripts", "make-data_woo2022_macrophage.R", package = "scpdata")`
##'
##' @source
##'
##' The peptide and protein data can be downloaded from the MASSIVE
##' repository MSV000085937
##'
##' @references
##'
##' **Source article**: Woo, Jongmin, Geremy C. Clair, Sarah M.
##' Williams, Song Feng, Chia-Feng Tsai, Ronald J. Moore, William B.
##' Chrisler, et al. 2022. “Three-Dimensional Feature Matching
##' Improves Coverage for Single-Cell Proteomics Based on Ion Mobility
##' Filtering.” Cell Systems 13 (5): 426–34.e4.
##' ([link to article](http://dx.doi.org/10.1016/j.cels.2022.02.003)).
##'
##' @examples
##' \donttest{
##' woo2022_macrophage()
##' }
##'
##' @keywords datasets
##'
"woo2022_macrophage"
