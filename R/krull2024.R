##' Krull et al, 2024 (Nature Communications): IFN-\eqn{\gamma} response
##'
##' They develop a new strategy for data-independent acquisition (DIA) that 
##' leverages the co-analysis of low-input samples alongside a corresponding 
##' enhancer (ME) of higher input. Using DIA-ME, they investigate the 
##' proteomic response of U-2 OS cells to interferon gamma (IFN-\eqn{\gamma}) at 
##' the single-cell level.
##'
##' @format A [QFeatures] object with 159 assays, each assay being a
##' [SingleCellExperiment] object.
##'
##' - Assay 1-158: DIA-NN main output report table split for each
##'   acquisition run. First 15 run acquires 10 single cells (MEs) and, 
##'   remaining 143 run acquires 1 single cell. It contains the results
##'   of the spectrum identification and quantification.
##' - `proteins`: DIA-NN protein group matrix, containing normalised
##'   quantities for 1553 protein groups in 143 single cells. Proteins
##'   are filtered at (Q.Value <= 0.01), (Lib.Q.Value <= 0.01), and 
##'   (Lib.PG.Q.Value <= 0.01).
##'   
##' The `colData(krull2024())` contains cell type annotations. The description 
##' of the `rowData` fields for the different assays can be found in the
##' [`DIA-NN` documentation](https://github.com/vdemichev/DiaNN#readme).
##'
##' @section Acquisition protocol:
##'
##' The data were acquired using the following setup. More information
##' can be found in the source article (see `References`).
##'
##' - **Cell isolation**: cells were detached with trypsin digestion, followed 
##'   by dilution in 1.5 mL PBS, and isolated using BD FACSAria III instrument.
##' - **Sample preparation**: Sorted single cells were collected in lysis 
##'   buffer (50 mM TEAB, pH 8.5, and 0.025% DDM), denatured at 70 degrees 
##'   Celsius for 30 minutes. Samples were acidified with 0.5% FA and 
##'   transferred to auto sampler plates for mass spectrometry analysis.
##' - **Separation**: Peptides were injected in a 2 microliter volume onto 
##'   a (25 cm x 75  micrometer) ID column at a flow rate of 300 nL/min, 
##'   separated using a gradient of ACN in water with 0.1% FA over 15 minutes, 
##'   connected to a nano-ESI source.
##' - **Ionization**: Ionization was performed using a 1,500 V capillary 
##'   voltage with 3.0 L/min dry gas and a dry temperature of 180 degrees 
##'   Celsius. MS data acquisition was conducted in diaPASEF mode using a 
##'   timsTOF Pro mass spectrometer.
##' - **Mass spectrometry**: MS1 scans covered a range of 200-1,700 m/z, 
##'   while DIA window isolation targeted 475-1,000 m/z with eight DIA scans 
##'   per cycle. Fragmentation was triggered by collision energy ranging from 
##'   45 eV to 27 eV depending on the ion mobility.
##' - **Data analysis**: Data was processed using DIA-NN (v1.8.0) and 
##'   Spectronaut 18 in a library-free approach, using deep learning 
##'   for spectrum prediction, retention times, and ion mobility.
##'
##' @section Data collection:
##'
##' The data were collected from the PRIDE
##' [repository](https://www.ebi.ac.uk/pride/archive/projects/PXD053464)
##' in the `03_SingleCell_Searches.zip` file.
##'
##' We loaded the DIA-NN main report table and generated a sample
##' annotation table based on the MS file names. We next combined the
##' sample annotation and the DIANN tables into a [QFeatures] object
##' following the `scp` data structure. We loaded the proteins group
##' matrix as a [SingleCellExperiment] object, and added the protein data 
##' as a new assay and link the precursors to proteins using the 
##' `Protein.Group` variable from the `rowData`.
##'
##' @source
##' The data were downloaded from PRIDE
##' [repository](https://www.ebi.ac.uk/pride/archive/projects/PXD053464)
##' with accession ID `PXD053464`.
##'
##' @references
##' Krull, K. K., Ali, S. A., & Krijgsveld, J. 2024. "Enhanced feature matching
##' in single-cell proteomics characterizes IFN-\eqn{\gamma} response and co-existence of
##' Cell States." Nature Communications, 15(1). 
##' [Link to article](https://doi.org/10.1038/s41467-024-52605-x) 
##' 
##' @examples
##' \donttest{
##' krull2024()
##' }
##'
##' @keywords datasets
##'
"krull2024"
