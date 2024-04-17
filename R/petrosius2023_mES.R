##' Petrosius et al, 2023 (Nat. Comm.): Mouse embryonic stem cell (mESC) in
##' different culture conditions
##'
##' @description
##' Profiling mouse embryonic stem cells across ground-state (m2i) and
##' differentiation-permissive (m15) culture conditions. The data were
##' acquired using orbitrap-based data-independent acquisition (DIA).
##' The objective was to demonstrate the capability of their approach
##' by profiling mouse embryonic stem cell culture conditions, showcasing
##' heterogeneity in global proteomes, and highlighting differences in
##' the expression of key metabolic enzymes in distinct cell subclusters.
##'
##' @format A [QFeatures] object with 605 assays, each assay being a
##' [SingleCellExperiment] object:
##'
##' - Assay 1-603: PSM data acquired with an orbitrap-based data-independent
##'   acquisition (DIA) protocol, hence those assays contain single column
##'   that contains the quantitative information.
##' - `peptides`: peptide data containing quantitative data for 9884
##'   peptides and 603 single-cells.
##' - `proteins`: protein data containing quantitative data for 4270
##'   proteins and 603 single-cells.
##'
##' Sample annotation is stored in `colData(petrosius2023_mES())`.
##'
##' @section Acquisition protocol:
##'
##' The data were acquired using the following setup. More information
##' can be found in the source article (see `References`).
##'
##' - **Sample isolation**: Cell sorting was done on a Sony MA900 cell sorter
##'   using a 130 microm sorting chip. Cells were sorted at single-cell resolution,
##'   into a 384-well Eppendorf LoBind PCR plate (Eppendorf AG) containing 1 microL
##'   of lysis buffer.
##' - **Sample preparation**: Single-cell protein lysates were digested with
##'   2 ng of Trypsin supplied in 1 microL of digestion buffer which was 
##'   carried out overnight at 37 °C, and subsequently acidified by the 
##'   addition of 1 microL 1% (v/v) trifluoroacetic acid (TFA). All liquid 
##'   dispensing was done using an I-DOT One instrument.
##' - **Liquid chromatography**: For the HRMS1-DIA experiments and the DIA 
##'   isolation window survey, the Evosep One liquid chromatography was used.
##'   The standard 31 min or 58min pre-defined Whisper gradients were used with a 
##'   flow rate of 100 nl/min for peptide elution. 
##' - **Mass spectrometry**: The mass spectrometer was operated in positive
##'   mode with the FAIMSPro interface compensation voltage set to -45 V.
##'   MS1 scans were carried out at 120,000 resolution with an automatic gain
##'   control (AGC) of 300% and maximum injection time set to auto. For the DIA
##'   isolation window survey a scan range of 500–900 was used and 400–1000
##'   rest of the experiments. Higher energy collisional dissociation (HCD) was
##'   used for precursor fragmentation with a normalized collision energy (NCE)
##'   of 33% and MS2 scan AGC target was set to 1000%.
##' - **Raw data processing**: The mESC raw data files were processed with
##'   Spectronaut 17.
##'
##' @section Data collection:
##'
##' The data were provided by the Author and is accessible at the
##' [Dataverse](https://dataverse.uclouvain.be/dataset.xhtml?persistentId=doi:10.14428/DVN/EMAVLT)
##' The folder ('20240205_111248_mESC_SNEcombine_m15-m2i/') contains
##' the following files of interest:
##'
##' - `20240205_111251_PEPQuant (Normal).tsv`: the PSM level data
##' - `20240205_111251_Peptide Quant (Normal).tsv`: the peptide level data
##' - `20240205_111251_PGQuant (Normal).tsv`: the protein level data
##'
##' The metadata were downloaded from the [Zenodo repository](https://zenodo.org/records/8146605).
##'
##' - `sample_facs.csv`: the metadata
##'
##' We formatted the quantification table so that columns match with the
##' metadata. Then, both tables are then combined in a single
##' [QFeatures] object using the [scp::readSCP()] function.
##'
##' The peptide data were formated to a [SingleCellExperiment] object and the
##' sample metadata were matched to the column names and stored in the `colData`.
##' The object is then added to the [QFeatures] object and the rows of the PSM
##' data are linked to the rows of the peptide data based on the peptide sequence
##' information through an `AssayLink` object.
##'
##' The protein data were formated to a [SingleCellExperiment] object and
##' the sample metadata were matched to the column names and stored in the
##' `colData`. The object is then added to the [QFeatures] object and the rows
##' of the peptide data are linked to the rows of the protein data based on the
##' protein sequence information through an `AssayLink` object.
##'
##' @source The peptide and protein data can be downloaded from the
##'     [Dataverse](https://dataverse.uclouvain.be/dataset.xhtml?persistentId=doi:10.14428/DVN/EMAVLT)
##'     The raw data and the quantification data can also be found in
##'     the MassIVE repository `MSV000092429`:
##'     ftp://MSV000092429@massive.ucsd.edu/.
##'
##' @references
##' **Source article**: Petrosius, V., Aragon-Fernandez, P., Üresin, N. et al.
##' "Exploration of cell state heterogeneity using single-cell proteomics
##' through sensitivity-tailored data-independent acquisition."
##' Nat Commun 14, 5910 (2023).
##' ([link to article](https://doi.org/10.1038/s41467-023-41602-1)).
##'
##' @examples
##' \donttest{
##' petrosius2023_mES()
##' }
##'
##' @keywords datasets
##'
"petrosius2023_mES"
