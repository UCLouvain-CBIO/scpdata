##' Zhu et al. 2018 (Nat. Comm.): human pancreatic islets
##'
##'
##' Near single-cell proteomics data human pancreas samples. The
##' samples were collected from pancreatic tissue slices using laser
##' dissection. The pancreata were obtained from organ donors through
##' the JDRFNetwork for Pancreatic Organ Donors with Diabetes (nPOD)
##' program. The sample come either from control patients (n=9) or
##' from type 1 diabetes (T1D) patients (n=9).
##'
##' @format A [QFeatures] object with 4 assays, each assay being a
##' [SingleCellExperiment] object:
##'
##' - `peptides`: quantitative information for 24,321 peptides from
##'   18 islet samples
##' - `proteins_intensity`: quantitative information for 3,278
##'   proteins from 18 islet samples
##' - `proteins_LFQ`: LFQ intensities for 3,278 proteins from 18 islet
##'   samples
##' - `proteins_iBAQ`: iBAQ values for 3,278 proteins from 18 islet
##'   samples
##'
##' Sample annotation is stored in `colData(zhu2018NC_islets())`.
##'
##' @section Acquisition protocol:
##'
##' The data were acquired using the following setup. More information
##' can be found in the source article (see `References`).
##'
##' - **Cell isolation**: The islets were extracted from the pacreatic
##'   tissues using laser-capture microdissection.
##' - **Sample preparation** performed using the nanoPOTs device.
##'   Protein extraction using RapiGest (+ DTT) + alkylation (IAA) +
##'   Lys-C digestion + cleave RapiGest (formic acid)
##' - **Separation**: nanoACQUITY UPLC pump with an Self-Pack PicoFrit
##'   70cm x 30um LC columns; 60nL/min)
##' - **Ionization**: ESI (1,900V)
##' - **Mass spectrometry**: Thermo Fisher Orbitrap Fusion Lumos
##'   Tribrid. MS1 settings: accumulation time = 246ms; resolution =
##'   120,000; AGC = 1E6. MS/MS settings: accumulation time = 118ms;
##'   resolution = 60,000; AGC = 1E5.
##' - **Data analysis**: MaxQuant (v1.5.3.30) + Perseus + OriginLab
##'   2017
##'
##' @section Data collection:
##'
##' The data were collected from the PRIDE repository (accession
##' ID: PXD006847).  We downloaded the `Islet_t1d_ct_peptides.txt`
##' and the `Islet_t1d_ct_proteinGroups.txt` files containing the
##' combined identification and quantification results. The sample
##' types were inferred from the names of columns holding the
##' quantification data. The peptides data were converted to a
##' [SingleCellExperiment] object. We split the protein table to
##' separate the three types of quantification: protein intensity,
##' label-free quantitification (LFQ) and intensity based absolute
##' quantification (iBAQ). Each table is converted to a
##' [SingleCellExperiment] object along with the remaining protein
##' annotations. The 4 objects are combined in a single [QFeatures]
##' object and feature links are created based on the peptide leading
##' razor protein ID and the protein ID.
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
##' @examples
##' \donttest{
##' zhu2018NC_islets()
##' }
##'
##' @keywords datasets
##'
"zhu2018NC_islets"
