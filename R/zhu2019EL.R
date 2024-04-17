##' Zhu et al. 2019 (eLife): chicken utricle cells
##'
##'
##' Single-cell proteomics data from chicken utricle acquired to
##' study the hair-cell development. The cells are isolated from
##' peeled utrical epithelium and separated into hair cells (FM1-43
##' high) and supporting cells (FM1-43 low). The sample contain either
##' 1 cell (n = 28), 3 cells (n = 7), 5 cells (n = 8) or 20 cells (n =
##' 14).
##'
##' @format A [QFeatures] object with 62 assays, each assay being a
##' [SingleCellExperiment] object:
##'
##' - `XYZw`: 60 assays containing PSM data. The sample are annotated
##'   as follows. `X` indicates the experiment, either 1 or 2. `Y`
##'   indicated the FM1-43 signal, either high (H) or low (L). `Z`
##'   indicates the number of cells (0, 1, 3, 5 or 20). `w` indicates
##'   the replicate, starting from `a`, it can go up to `j`.
##' - `peptides`: quantitative data for 3444 peptides in 60 samples
##'   (all runs are combined).
##' - `proteins_intensity`: protein intensities for 840 proteins
##'   from 24 samples
##' - `proteins_iBAQ`: iBAQ values for 840 proteins from 24 samples
##'
##' Sample annotation is stored in `colData(zhu2019EL())`.
##'
##' @section Acquisition protocol:
##'
##' The data were acquired using the following setup. More information
##' can be found in the source article (see `References`).
##'
##' - **Cell isolation**: The cells were taken from the utricles of
##'   E15 chick embryos. Samples were stained with FM1-43FX and the
##'   cells were dissociated using enzymatic digestion. Cells were
##'   FACS sorted (BD Influx) and split based on their FM1-43 signal,
##'   while ensuring no debris, doublets or dead cells are retained.
##' - **Sample preparation** performed using the nanoPOTs device. Cell
##'   lysis and protein extraction and reduction are performed using
##'   dodecyl beta-D-maltoside + DTT + ammonium bicarbonate. Protein
##'   were then alkylated using IAA. Protein digestion is performed
##'   using Lys-C and trypsin. Finally samples acidification is
##'   performed using formic acid.
##' - **Separation**: Dionex UltiMate pump with an C18-Packed column
##'   (50cm x 30um; 60nL/min)
##' - **Ionization**: ESI (2,000V)
##' - **Mass spectrometry**: Orbitrap Fusion Lumos Tribrid. MS1
##'   settings: accumulation time = 246ms; resolution = 120,000; AGC =
##'   3E6. MS/MS settings: accumulation time = 502ms; resolution =
##'   120,000; AGC = 2E5.
##' - **Data analysis**: Andromeda & MaxQuant (v1.5.3.30) and the
##'   search database is NCBI GRCg6a.
##'
##' @section Data collection:
##'
##' All data were collected from the PRIDE repository (accession ID:
##' PXD014256).
##'
##' The sample annotation information is provided in the
##' `Zhu_2019_chick_single_cell_samples_CORRECTED.xlsx` file. This file
##' was given during a personal discussion and is a corrected version
##' of the annotation table available on the PRIDE repository.
##'
##' The PSM data were found in the `evidence.txt` (in the
##' `Experiment 1+ 2`) folder. The PSM data were filtered so that it
##' contains only samples that are annotated. The data were then
##' converted to a [QFeatures] object using the [scp::readSCP()]
##' function.
##'
##' The peptide data were found in the `peptides.txt` file. The column
##' names holding the quantitative data were adapted to match the
##' sample names in the [QFeatures] object. The data were then
##' converted to a [SingleCellExperiment] object and then inserted in
##' the [QFeatures] object. Links between the PSMs and the peptides
##' were added
##'
##' A similar procedure was applied to the protein data. The data were
##' found in the `proteinGroups.txt` file. We split the protein table
##' to separate the two types of quantification: summed intensity and
##' intensity based absolute quantification (iBAQ). Both tables are
##' converted to [SingleCellExperiment] objects and are added to the
##' [QFeatures] object as well as the `AssayLink` between peptides and
##' proteins.
##'
##' @source
##' The PSM data can be downloaded from the PRIDE repository
##' PXD014256. The source link is:
##' ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/11/PXD014256
##'
##' @references
##'
##' Zhu, Ying, Mirko Scheibinger, Daniel Christian Ellwanger, Jocelyn
##' F. Krey, Dongseok Choi, Ryan T. Kelly, Stefan Heller, and Peter G.
##' Barr-Gillespie. 2019. “Single-Cell Proteomics Reveals Changes in
##' Expression during Hair-Cell Development.” eLife 8 (November).
##' ([link to article](https://doi.org/10.7554/eLife.50777)).
##'
##' @examples
##' \donttest{
##' zhu2019EL()
##' }
##'
##' @keywords datasets
##'
"zhu2019EL"
