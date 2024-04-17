##' Brunner et al. 2022 (Mol. Syst. Biol.): cell cycle state study
##'
##' Single cell proteomics data acquired by the Mann Lab using a newly
##' designed timsTOF instrument, referred to as timsTOF-SCP. The
##' dataset contains quantitative information from single-cells blocked
##' at 4 cell cycle stages: G1, G1-S, G2, G2-M. The data was acquired
##' using a label-free sample preparation protocole combined to a
##' data independent (DIA) acquisition mode.
##'
##' @format A [QFeatures] object with 435 assays, each assay being a
##' [SingleCellExperiment] object.
##'
##' - Assay 1-434: DIA-NN main output report table split for each
##'   acquisition run. Since each run acquires 1 single cell, each
##'   assay contains a single column. It contains the results
##'   of the spectrum identification and quantification.
##' - `protein`: DIA-NN protein group matrix, containing normalised
##'   quantities for 2476 protein groups in 434 single cells. Proteins
##'   are filtered at 1% FDR, using global q-values for protein groups
##'   and both global and run-specific q-values for precursors.
##'
##' The `colData(brunner2022())` contains cell type annotations and
##' batch annotations. The description of the `rowData` fields for the
##' different assays can be found in the
##' [`DIA-NN` documentation](https://github.com/vdemichev/DiaNN#readme).
##'
##' @section Acquisition protocol:
##'
##' The data were acquired using the following setup. More information
##' can be found in the source article (see `References`).
##'
##' - **Cell isolation**: cells were detached with trypsin treatment,
##'   followed by strong pipetting, and isolate using FACS.
##' - **Sample preparation**: cell lysis by freeze-heat followed by
##'   sonication, overnight protein digestion with trypsin/lysC mix and
##'   desalting using EvoTips trap column (EvoSep)
##' - **Separation**: online EvoSep One LC system using a 5 cm x 75 µm
##'   ID column with 1.9µm C18 beads (EvoSep) at 100nL/min flow rate.
##' - **Ionization**: 10µm ID zero dead volume electrospray emitter
##'   (Bruker Daltonik) + nanoelectro-spray ion source (Captive spray,
##'   Bruker Daltonik)
##' - **Mass spectrometry**: DIA PASEF mode. Correlation between IM
##'   and m/z was used to synchronize the elution of precursors from
##'   each IM scan with the quadrupole isolation window. Five
##'   consecutive diaPASEF cycles. The collision energy was ramped
##'   linearly as a function of the IM from 59 eV at 1/K0=1.6 Vs cm^2
##'   to 20 eV at 1/K0=0.6 Vs cm^2.
##' - **Data analysis**: DIA-NN (1.8).
##'
##' @section Data collection:
##'
##' The data were collected from the PRIDE
##' [repository](https://www.ebi.ac.uk/pride/archive/projects/PXD024043)
##' in the `DIANN1.8_SingleCells_CellCycle.zip` file.
##'
##' We loaded the DIA-NN main report table and generated a sample
##' annotation table based on the MS file names. We next combined the
##' sample annotation and the DIANN tables into a [QFeatures] object
##' following the `scp` data structure. We loaded the proteins group
##' matrix as a [SingleCellExperiment] object, fixed ambiguous
##' protein group names, and added the protein data as a new assay and
##' link the precursors to proteins using the `Protein.Group` variable
##' from the `rowData`.
##'
##' @source
##' The data were downloaded from PRIDE
##' [repository](https://www.ebi.ac.uk/pride/archive/projects/PXD024043)
##' with accession ID `PXD024043`.
##'
##' @references
##' Brunner, Andreas-David, Marvin Thielert, Catherine Vasilopoulou,
##' Constantin Ammar, Fabian Coscia, Andreas Mund, Ole B. Hoerning, et
##' al. 2022. "Ultra-High Sensitivity Mass Spectrometry Quantifies
##' Single-Cell Proteome Changes upon Perturbation." Molecular Systems
##' Biology 18 (3): e10798.
##' [Link to article](http://dx.doi.org/10.15252/msb.202110798)
##'
##' @examples
##' \donttest{
##' brunner2022()
##' }
##'
##' @keywords datasets
##'
"brunner2022"
