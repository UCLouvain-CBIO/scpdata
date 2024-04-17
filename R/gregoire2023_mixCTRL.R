##' Grégoire et al. 2023 - mixCTRL (arXiv): benchmark using
##' monocytes/macrophages
##'
##' Single cell proteomics data acquired using the SCoPE2 protocol.
##' The dataset contains two monocytes cell lines (THP1 and U937) as
##' well as controled mixtures of both and macrophage-like cells
##' produced upon PMA treatment. It contains quantitative information
##' at PSM, peptide and protein levels. Data was acquired using Lumos
##' Orbitrap (mainly) and timsTOF SCP mass spectrometers.
##'
##' @format A [QFeatures] object with 119 assays, each assay being a
##' [SingleCellExperiment] object:
##'
##' - Assays 1-42: PSM data acquired with a TMT-16plex protocol, hence
##'   those assays contain 16 columns. Columns hold quantitative
##'   information from single-cell channels, carrier channels,
##'   blank (negative control) channels and unused channels.
##' - Assays 43-84: peptide data resulting from the PSM to peptide
##'   aggregation of the 42 PSM assays.
##' - Assays 85-91: peptide data for each of the 7 acquisition
##'   batches. Peptide data were joined based on their respective
##'   acquisition batches.
##' - Assays 92-98: normalised peptide data.
##' - Assays 99-105: normalised and log-transformed peptide data.
##' - Assays 106-112: protein data for each of the 7 acquisition
##'   batches. Normalised and log-transformed peptide data were
##'   agreggated to protein.
##' - Assays 113-119: Batch corrected protein data. Normalised and
##'   log-transformed protein data were batch corrected to remove
##'   technical variability induced by runs and channels.
##'
##' All the data has been filtered to keep high quality features and
##' samples.
##'
##' The `colData(gregoire2023_mixCTRL())` contains cell type annotation and
##' batch annotation that are common to all assays. The description of
##' the `rowData` fields for the PSM data can be found in the
##' [`sage` documentation](https://sage-docs.vercel.app/docs/results/search).
##'
##' @section Acquisition protocol:
##'
##' The data were acquired using the following setup. More information
##' can be found in the source article (see *References*).
##'
##' - **Cell isolation**: BD FACSAria III cell sorting.
##' - **Sample preparation** performed using the SCoPE2 protocol: mPOP
##'   cell lysis + trypsin digestion + TMT-16plex labeling and
##'   pooling.
##' - **Separation**: online nLC (Ultimate 3000 LC System or Vanquish
##'   Neo UHPLC System) with a BioZen Peptide Polar C18 250 x 0.0075mm
##'   column.
##' - **Mass spectrometry**:  Orbitrap Fusion Lumos Tribrid (MS1
##'   resolution = 70,000; MS2 accumulation time = 120ms; MS2
##'   resolution = 70,000) and timsTOF SCP.
##' - **Data preprocessing**: Sage.
##'
##' @section Data collection:
##'
##' The PSM data were collected from a Zenodo archive (see `Source`
##' section). The folder contains the following files of interest:
##'
##' - `results.sage.cbio.tsv`: the sage identification output file for
##'   batches acquired on the Lumos MS.
##' - `results.sage.giga.tsv`: the sage identification output file for
##'   batches acquired on the timsTOF SCP MS.
##' - `quant.cbio.tsv`: the sage quantification output file for
##'   batches acquired on the Lumos MS.
##' - `quant.giga.tsv`: the sage quantification output file for
##'   batches acquired on the timsTOF SCP MS.
##' - `sampleAnnotation_batch.csv`: sample annotation for each
##'   acquisition batch. There are in total 8 different annotation
##'   files.
##'
##' We combined the sample annotations in a single table. We also
##' combined `cbio` and `giga` tables together and merged resulting
##' identification and quantification tables. Both annotation and
##' features tables are then combined in a single [QFeatures] object
##' using the [scp::readSCP()] function.
##'
##' The [QFeatures] object was processed as described in the author's
##' manuscript (see `source`). Note that the imputed assays were used
##' in the paper for illustrative purposes only and have not been
##' reproduced here.
##'
##' @source
##' The data were downloaded from the [Zenodo
##' repository](https://zenodo.org/records/8417228).  The raw data and
##' the quantification data can also be found in the ProteomeXchange
##' Consortium via the [PRIDE partner
##' repository](https://www.ebi.ac.uk/pride/archive/projects/PXD046211),
##' project `PXD046211`.
##'
##' @references
##' Samuel Grégoire, Christophe Vanderaa, Sébastien Pyr dit Ruys,
##' Gabriel Mazzucchelli, Christopher Kune, Didier Vertommen and
##' Laurent Gatto. 2023. *Standardised workflow for mass spectrometry-
##' based single-cell proteomics data processing and analysis using
##' the scp package.*
##' arXiv. DOI:[10.48550/arXiv.2310.13598](https://doi.org/10.48550/arXiv.2310.13598)
##'
##' @examples
##' \donttest{
##' gregoire2023_mixCTRL()
##' }
##'
##' @keywords datasets
##'
"gregoire2023_mixCTRL"
