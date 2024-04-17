##' Guise et al. 2020 (Cell Rep.): postmortem ALS spinal moto neurons
##'
##' Single-cell proteomics data from postmortem human spinal moto
##' neurons (MN) obtained from control donors or donors with amyotrophic
##' lateral sclerosis (ALS). The data were generated following the
##' NanoPOTS protocol. Cells were isolated from samples obtained by
##' the university of Miami Brain Bank using laser capture
##' microdissection (LCM). Additional information about the amount of
##' TDP-43 intra-cellular levels has been assigned into levels 0 to 4.
##'
##' @format A [QFeatures] object with 102 assays, each assay being a
##' [SingleCellExperiment] object:
##'
##' - `F*`: 100 assays containing PSM data.
##' - `peptides`: quantitative data for 34,315 peptides in 108 samples.
##'   All samples combined, along with 8 additional unannotated
##'   samples.
##' - `proteins`: quantitative data for 4,437 protein groups in 108
##'   samples. All samples combined, along with 8 additional
##'   unannotated samples.
##'
##' Sample annotation is stored in `colData(guise2024())`.
##'
##' @section Acquisition protocol:
##'
##' The data were acquired using the following setup. More information
##' can be found in the source article (see `References`).
##'
##' - **Cell isolation**: The MN were isolated from samples obtained
##'   by the university of Miami Brain Bank using LCM.
##' - **Sample preparation** performed using the nanoPOTS workflow.
##'   Cells are treated with 0.1% DDM (for lysis) added with DTT
##'   (protein reduction), then IAA (alkylation), then Lys-C and
##'   trypsin (protein digestion).
##' - **Separation**: Samples were injected on the column using an
##'   Ultimate 3000 RSLCnano pump. The in-line loading column is a
##'   home-packed SPE column (5cm x 75um) while the peptide
##'   separation is performed on a an in-house-packed analytical SPE
##'   column (50 cm x 30um), using a 20nL/min flow rate.
##' - **Ionization**: nanospray emmitter (2,000V)
##' - **Mass spectrometry**: Orbitrap Exploris 480. HCD fragmentation.
##'   MS1 settings: accumulation time = 200 ms; resolution = 120,000;
##'   AGC = 1E6. MS2 settings: exclusion duration = 90 s;
##'   accumulation time = 500 ms; resolution = 30,000; AGC = 1E5.
##' - **Data analysis**: Sequest HT in Proteome Discoverer (v2.5) and
##'   the search database is Swiss-Prot (July 2020).
##'
##' @section Data collection:
##'
##' All data were collected from the MassIVE repository (accession ID:
##' MSV000092119).
##'
##' The sample annotations were combined from the tables in
##' `Biogen_TDP43_Round2_Reanalysis_10-13-2021_InputFiles.txt` and in
##' `Groups.txt`.
##'
##' The PSM data were found in the
##' `Biogen_TDP43_Round2_Reanalysis_10-13-2021_PSMs.txt` file. The
##' data were converted to a [QFeatures] object using the [scp::readSCP()]
##' function. We could not find sample annotations for MS run ID:
##' F61, F34, F42, F88, F77, F8, F21, F5.
##'
##' The peptide data were found in the
##' `Biogen_TDP43_Round2_Reanalysis_10-13-2021_PeptideGroups.txt`
##' file. The column names holding the quantitative data were adapted
##' to match the sample names in the [QFeatures] object. The data were
##' then converted to a [SingleCellExperiment] object and then
##' inserted in the [QFeatures] object.
##'
##' A similar procedure was applied to the protein data. The data were
##' found in the
##' `Biogen_TDP43_Round2_Reanalysis_10-13-2021_Proteins.txt` file. The
##'  column names were
##' adapted, the data were converted to a [SingleCellExperiment]
##' object and then inserted in the [QFeatures] object.
##'
##' @source
##'
##' All data can be downloaded from the MassIVE repository
##' MSV000092119. The source link is:
##' ftp://massive.ucsd.edu/v05/MSV000092119/
##'
##' @references
##'
##' Guise, Amanda J., Santosh A. Misal, Richard Carson, Jen-Hwa Chu,
##' Hannah Boekweg, Daisha Van Der Watt, Nora C. Welsh, et al. 2024.
##' “TDP-43-Stratified Single-Cell Proteomics of Postmortem Human
##' Spinal Motor Neurons Reveals Protein Dynamics in Amyotrophic
##' Lateral Sclerosis.” Cell Reports 43 (1): 113636.
##' ([link to article](http://dx.doi.org/10.1016/j.celrep.2023.113636)).
##'
##' @examples
##' \donttest{
##' guise2024()
##' }
##'
##' @keywords datasets
##'
"guise2024"
