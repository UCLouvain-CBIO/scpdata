##' Zhu et al. 2018 (Nat. Comm.): HeLa lysates
##'
##' Near single-cell proteomics data of HeLa lysates at different
##' concentrations (10, 40 and 140 cell equivalent). Each
##' concentration is acquired in triplicate.
##'
##' @format A [QFeatures] object with 4 assays, each assay being a
##' [SingleCellExperiment] object:
##'
##' - `peptides`: quantitative information for 14,921 peptides from
##'   9 lysate samples
##' - `proteins_intensity`: quantitative information for 2,199
##'   proteins from 9 lysate samples
##' - `proteins_LFQ`: LFQ intensities for 2,199 proteins from 9 lysate
##'   samples
##' - `proteins_iBAQ`: iBAQ values for 2,199 proteins from 9 lysate
##'   samples
##'
##' Sample annotation is stored in `colData(zhu2018NC_lysates())`.
##'
##' @section Acquisition protocol:
##'
##' The data were acquired using the following setup. More information
##' can be found in the original article (see `References`).
##'
##' - **Cell isolation**: HeLas were collected from cell cultures.
##' - **Sample preparation** performed in bulk (5E5 cells/mL). Protein
##'   extraction using RapiGest (+ DTT) + dilution to target
##'   concentration + alkylation (IAA) + Lys-C digestion + trypsin
##'   digestion + cleave RapiGest (formic acid).
##' - **Separation**: nanoACQUITY UPLC pump (60nL/min) with an
##'   Self-Pack PicoFrit 70cm x 30um LC columns.
##' - **Ionization**: ESI (1,900V).
##' - **Mass spectrometry**: Thermo Fisher Orbitrap Fusion Lumos
##'   Tribrid. MS1 settings: accumulation time = 246ms; resolution =
##'   120,000; AGC = 1E6. MS/MS settings, depend on the sample size,
##'   excepted for the AGC = 1E5. Blank and approx. 10 cells (time = 502ms;
##'   resolution = 240,000), approx. 40 cells (time = 246ms; resolution =
##'   120,000), approx. 140 cells (time = 118ms; resolution = 60,000).
##' - **Data analysis**: MaxQuant (v1.5.3.30) + Perseus + OriginLab
##'   2017.
##'
##' @section Data collection:
##'
##' The data were collected from the PRIDE repository (accession
##' ID: PXD006847).  We downloaded the `Vail_Prep_Vail_peptides.txt`
##' and the `Vail_Prep_Vail_proteinGroups.txt` files containing the
##' combined identification and quantification
##' results. The sample annotations were inferred from the names of
##' columns holding the quantification data and the information in the
##' article. The peptides data were converted to a [SingleCellExperiment]
##' object. We split the protein table to separate the three types of
##' quantification: protein intensity, label-free quantitification
##' (LFQ) and intensity based absolute quantification (iBAQ). Each
##' table is converted to a [SingleCellExperiment] object along with
##' the remaining protein annotations. The 4 objects are combined in
##' a single [QFeatures] object and feature links are created based on
##' the peptide leading razor protein ID and the protein ID.
##'
##' @source
##' The PSM data can be downloaded from the PRIDE repository
##' PXD006847. The source link is:
##' ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2018/01/PXD006847
##'
##' @references
##'
##' Zhu, Ying, Paul D. Piehowski, Rui Zhao, Jing Chen, Yufeng Shen,
##' Ronald J. Moore, Anil K. Shukla, et al. 2018. “Nanodroplet
##' Processing Platform for Deep and Quantitative Proteome Profiling
##' of 10-100 Mammalian Cells.” Nature Communications 9 (1): 882
##' ([link to article](http://dx.doi.org/10.1038/s41467-018-03367-w)).
##'
##' @seealso The same experiment was conducted directly on HeLa cells
##' samples rather than lysates. The data is available in
##' [zhu2018NC_hela].
##'
##' @examples
##' \donttest{
##' zhu2018NC_lysates()
##' }
##'
##' @keywords datasets
##'
##'
"zhu2018NC_lysates"
