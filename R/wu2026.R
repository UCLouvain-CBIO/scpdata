##' Wu et al., 2026 (Nature Biotechnology): single-cell proteomic landscape
##' of the developing human brain
##'
##' @description
##'
##' Label-free single-cell proteomics data from prenatal human brain tissue
##' (gestational weeks 13, 15 and 19), generated to characterize cell-type
##' proteomes and developmental trajectories at single-cell resolution.
##'
##' @format A [QFeatures] object with 2313 sets, each set being a
##' [SingleCellExperiment] object:
##'
##' - Set 1-2310: log-transformed precursor-level single-cell data (one
##'   column per set) from DIA-NN reports.
##' - `peptides`: peptide-level data containing quantitative data for 17755
##'   peptides in 2310 single cells.
##' - `imported_proteins`: protein-level table imported from the authors'
##'   protein report, containing 4401 proteins in 2310 single cells.
##' - `proteins`: protein-level data containing quantitative data for 5107
##'   proteins in 2310 single cells.
##'
##' Sample annotation is stored in `colData(wu2026())`.
##'
##' @section Acquisition protocol:
##'
##' The data were acquired using the following setup. More information
##' can be found in the source article (see `References`).
##'
##' - **Sample isolation**: fresh prenatal human brain tissue was dissociated
##'   and single cells were sorted by FACS (Sony SH800) into 384-well low-bind
##'   plates.
##' - **Sample preparation**: cells were lysed in low-volume TEAB-based buffer,
##'   digested overnight with trypsin, acidified with trifluoroacetic acid,
##'   and desalted on EvoTip Pure.
##' - **Liquid chromatography**: peptides were separated on an IonOpticks Aurora
##'   Elite XT C18 nano-LC column (15 cm x 75 micrometer) on an Evosep One
##'   using a 31-minute gradient.
##' - **Mass spectrometry**: Orbitrap Eclipse Tribrid with FAIMS Pro was used
##'   in DIA mode.
##' - **Raw data processing**: DIA-NN against the UniProt human reference
##'   proteome.
##'
##' @section Data collection:
##'
##' The `QFeatures` object was built from the DIA-NN peptide and protein report
##' tables and a sample annotation table. PSM-level data were imported with
##' [scp::readSCP()], log-transformed, aggregated to peptide-level features, and
##' joined across single-cell runs. The protein report table was imported as a
##' [SingleCellExperiment], added as an set, and linked to peptide-level data.
##'
##' The script to reproduce the `QFeatures` object is available at
##' `inst/scripts/make-data_wu2026.R`.
##'
##' @source
##' The data were downloaded from the PRIDE repository with accession ID
##' `PXD071075`:
##' https://www.ebi.ac.uk/pride/archive/projects/PXD071075
##'
##' @references
##' Wu, T., Jiang, L., Mukhtar, T., Wang, L., Jian, R., Wang, C., Trinh, T.,
##' Kriegstein, A. R., Snyder, M., and Li, J. 2026. "Single-cell proteomic
##' landscape of the developing human brain." Nature Biotechnology.
##' [Link to article](https://doi.org/10.1038/s41587-025-02980-7)
##'
##' @examples
##' \donttest{
##' wu2026()
##' }
##'
##' @keywords datasets
##'
"wu2026"
