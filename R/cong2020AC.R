##' Cong et al. 2020 (Ana. Chem.): HeLa single cells
##'
##' Single-cell proteomics using the nanoPOTS sample processing device
##' in combination with ultranarrow-bore (20um i.d.) packed-column LC
##' separations and the Orbitrap Eclipse Tribrid MS. The dataset
##' contains label-free quantitative information at PSM, peptide and
##' protein level. The samples are single Hela cells. Bulk samples
##' (100 and 20 cells) were also included in the experiment to
##' increase the idendtification rate thanks to between-run matching
##' (cf MaxQuant).
##'
##' @format A [QFeatures] object with 9 assays, each assay being a
##' [SingleCellExperiment] object:
##'
##' - `100/20 HeLa cells`: 2 assays containing PSM data for a bulk
##'   of 100 or 20 HeLa cells, respectively.
##' - `Blank`: assay containing the PSM data for a blank sample
##' - `Single cell X`: 4 assays containing PSM data for a single cell.
##'   The `X` indicates the replicate number.
##' - `peptides`: quantitative data for 12590 peptides in 7 samples
##'   (all runs combined).
##' - `proteins`: quantitative data for 1801 proteins in 7 samples
##'   (all runs combined).
##'
##' Sample annotation is stored in `colData(cong2020AC())`.
##'
##' @section Acquisition protocol:
##'
##' The data were acquired using the following setup. More information
##' can be found in the source article (see `References`).
##'
##' - **Cell isolation**: The HeLa cells were diluted and aspired
##'   using a microcapillary with a pulled tip.
##' - **Sample preparation** performed using the nanoPOTs device.
##'   Protein extraction using RapiGest (+ DTT) + alkylation (IAA) +
##'   Lys-C digestion + cleave RapiGest (formic acid)
##' - **Separation**: UltiMate 3000 RSLCnano pump with a home-packed
##'   nanoLC column (60cm x 20um i.d.; approx. 20 nL/min)
##' - **Ionization**: ESI (2,000V; Nanospray Flex)
##' - **Mass spectrometry**: Thermo Fisher Orbitrap Fusion Eclipse.
##'   MS1 settings: accumulation time = 246ms; resolution = 120,000;
##'   AGC = 1E6. MS/MS settings depend on quantity. All: AGC = 1E5.
##'   20-100 cels: accumulation time = 246ms; resolution = 120,000.
##'   Single cells: accumulation time = 500ms; resolution = 240,000.
##' - **Data analysis**: MaxQuant (v1.6.3.3) + Excel
##'
##' @section Data collection:
##'
##' The PSM, peptide and protein data were collected from the PRIDE
##' repository (accession ID: PXD016921).  We downloaded the
##' `evidence.txt` file containing the PSM identification and
##' quantification results. The sample annotation was inferred from
##' the samples names. The data were then converted to a [QFeatures]
##' object using the [scp::readSCP()] function.
##'
##' The peptide data were processed similarly from the `peptides.txt`
##' file. The quantitative column names were adpated to match the PSM
##' data. The peptide data were added to [QFeatures] object and link
##' between the features were stored.
##'
##' The protein data were similarly processed from the
##' `proteinGroups.txt` file. The quantitative column names were
##' adapted to match the PSM data. The peptide data were added to
##' [QFeatures] object and link between the features were stored.
##'
##' @source
##' All files can be downloaded from the PRIDE repository PXD016921.
##' The source link is:
##' ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2020/02/PXD016921
##'
##' @references
##'
##' Cong, Yongzheng, Yiran Liang, Khatereh Motamedchaboki, Romain
##' Huguet, Thy Truong, Rui Zhao, Yufeng Shen, Daniel Lopez-Ferrer,
##' Ying Zhu, and Ryan T. Kelly. 2020. “Improved Single-Cell Proteome
##' Coverage Using Narrow-Bore Packed NanoLC Columns and
##' Ultrasensitive Mass Spectrometry.” Analytical Chemistry, January.
##' ([link to article](https://doi.org/10.1021/acs.analchem.9b04631)).
##'
##' @examples
##' \donttest{
##' cong2020AC()
##' }
##'
##' @keywords datasets
##'
"cong2020AC"
