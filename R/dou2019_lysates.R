##' Dou et al. 2019 (Anal. Chem.): HeLa lysates
##'
##' @description
##'
##' Single-cell proteomics using nanoPOTS combined with TMT
##' multiplexing. It contains quantitative information at PSM and
##' protein level. The samples are commercial Hela lysates diluted to
##' single-cell amounts (0.2 ng). The boosting wells contain the same
##' digest but at higher amount (10 ng).
##'
##' @format A [QFeatures] object with 3 assays, each assay being a
##' [SingleCellExperiment] object:
##'
##' - `Hela_run_1`: PSM data with 10 columns corresponding to the
##'   TMT-10plex channels. Columns hold quantitative information for
##'   HeLa lysate samples (either 0, 0.2 or 10ng). This is the data
##'   for run 1.
##' - `Hela_run_1`: PSM data with 10 columns corresponding to the
##'   TMT-10plex channels. Columns hold quantitative information for
##'   HeLa lysate samples (either 0, 0.2 or 10ng). This is the data
##'   for run 2.
##' - `peptides`: peptide data containing quantitative data for 13,934
##'   peptides in 20 samples (run 1 and run 2 combined).
##' - `proteins`: protein data containing quantitative data for 1641
##'   proteins in 20 samples (run 1 and run 2 combined).
##'
##' Sample annotation is stored in `colData(dou2019_lysates())`.
##'
##' @section Acquisition protocol:
##'
##' The data were acquired using the following setup. More information
##' can be found in the source article (see `References`).
##'
##' - **Cell isolation**: commercially available HeLa protein digest
##'   (Thermo Scientific).
##' - **Sample preparation** performed using the nanoPOTs device.
##'   Protein extraction (DMM + TCEAP) + alkylation (IAA) + Lys-C
##'   digestion + trypsin digestion + TMT-10plex labeling and pooling.
##' - **Separation**: nanoLC (Dionex UltiMate with an in-house packed
##'   50cm x 30um LC columns; 50nL/min)
##' - **Ionization**: ESI (2,000V)
##' - **Mass spectrometry**: Thermo Fisher Orbitrap Fusion Lumos
##'   Tribrid (MS1 accumulation time = 50ms; MS1 resolution = 120,000;
##'   MS1 AGC = 1E6; MS2 accumulation time = 246ms; MS2 resolution =
##'   60,000; MS2 AGC = 1E5)
##' - **Data analysis**: MS-GF+ + MASIC (v3.0.7111) + RomicsProcessor
##'   (custom R package)
##'
##' @section Data collection:
##'
##' The PSM data were collected from the MassIVE repository
##' MSV000084110 (see `Source` section). The downloaded files are:
##'
##' - `Hela_run_*_msgfplus.mzid`: the MS-GF+ identification result
##'   files
##' - `Hela_run_*_ReporterIons.txt`: the MASIC quantification result
##'   files
##'
##' For each batch, the quantification and identification data were
##' combined based on the scan number (common to both data sets). The
##' combined datasets for the different runs were then concatenated
##' feature-wise. To avoid data duplication due to ambiguous matching
##' of spectra to peptides or ambiguous mapping of peptides to proteins,
##' we combined ambiguous peptides to peptides groups and proteins to
##' protein groups. Feature annotations that are not common within a
##' peptide or protein group are are separated by a `;`. The sample
##' annotation table was manually created based on the available
##' information provided in the article. The data were then converted
##' to a [QFeatures] object using the [scp::readSCP()] function.
##'
##' We generated the peptide data. First, we removed PSM matched to
##' contaminants or decoy peptides and ensured a 1% FDR. We aggregated
##' the PSM to peptides based on the peptide (or peptide group)
##' sequence(s) using the median PSM instenity. The peptide data for
##' the different runs were then joined in a single assay (see
##' [QFeatures::joinAssays]), again based on the peptide sequence(s).
##' We then removed the peptide groups. Links between the peptide and
##' the PSM data were created using [QFeatures::addAssayLink]. Note
##' that links between PSM and peptide groups are not stored.
##'
##' The protein data were downloaded from `Supporting information`
##' section from the publisher's website (see `Sources`). The data is
##' supplied as an Excel file `ac9b03349_si_003.xlsx`. The file
##' contains 7 sheets from which we only took the sheet 6 (named
##' `5 - Run 1 and 2 raw data`) with the combined protein data for the
##' two runs. We converted the data to a [SingleCellExperiment]
##' object and added the object as a new assay in the [QFeatures]
##' dataset (containing the PSM data). Links between the proteins and
##' the peptides were created. Note that links to protein groups are
##' not stored.
##'
##' @source
##' The PSM data can be downloaded from the massIVE repository
##' MSV000084110. FTP link: ftp://massive.ucsd.edu/MSV000084110/
##'
##' The protein data can be downloaded from the
##' [ACS Publications](https://pubs.acs.org/doi/10.1021/acs.analchem.9b03349)
##' website (Supporting information section).
##'
##' @references
##' Dou, Maowei, Geremy Clair, Chia-Feng Tsai, Kerui Xu, William B.
##' Chrisler, Ryan L. Sontag, Rui Zhao, et al. 2019. “High-Throughput
##' Single Cell Proteomics Enabled by Multiplex Isobaric Labeling in a
##' Nanodroplet Sample Preparation Platform.” Analytical Chemistry,
##' September
##' ([link to article](https://doi.org/10.1021/acs.analchem.9b03349)).
##'
##' @seealso
##' [dou2019_mouse], [dou2019_boosting]
##'
##' @examples
##' \donttest{
##' dou2019_lysates()
##' }
##'
##' @keywords datasets
##'
##'
"dou2019_lysates"
