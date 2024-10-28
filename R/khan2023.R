##' Khan et al, 2023 (biorRxiv): Epithelial–Mesenchymal Transition
##'
##' @description
##'
##' Single-cell samples were prepared using the nPOP sample
##' preparation method.  Proteomics data were acquired using the
##' SCoPE2 protocol on a Thermo Scientific Q-Exactive mass
##' spectrometer. The dataset contains quantitative information on 421
##' MCF-10A single cells undergoing epithelial–mesenchymal transition
##' (EMT) triggered by TGF beta. The data are available at the PSM,
##' and protein levels. The paper investigates the dynamics of
##' correlation modules at the protein level.
##'
##' @format A [QFeatures] object with 47 assays, each assay being a
##' [SingleCellExperiment] object:
##'
##' - Assay 1-44: PSM data acquired with a TMTPro 16plex protocol, hence
##'   those assays contain 16 columns. Columns hold quantitative information
##'   from single-cell channels, carrier channels, reference channels,
##'   empty (negative control) channels and unused channels.
##' - `peptides`: peptide data containing quantitative data for 10055
##'   peptides and 421 single-cells.
##' - `proteins_imputed`: protein data containing quantitative data for 4571
##'   proteins and 420 single-cells with k-nearest neighbors (KNN) imputation.
##' - `proteins_unimputed`: protein data containing quantitative data for 4571
##'   proteins and 420 single-cells without imputation.
##'
##' The `colData(khan2023())` contains cell type and batch annotations that
##' are common to all assays. The description of the `rowData` fields for the
##' PSM data can be found in the
##' [`MaxQuant` documentation](https://cox-labs.github.io/coxdocs/output_tables.html).
##'
##' @section Acquisition protocol:
##'
##' The data were acquired using the following setup. More information
##' can be found in the source article (see `References`).
##'
##' - **Cell isolation**: CellenONE cell sorting.
##' - **Sample preparation** performed using the SCoPE2 protocol. nPOP
##'   cell lysis (DMSO) + trypsin digestion + TMTPro 16plex protocol.
##' - **Separation**: online nLC (DionexUltiMate 3000 UHPLC with a
##'   25cm x 75um IonOpticks Odyssey Series column (ODY3-25075C18); 200nL/min).
##' - **Ionization**: ESI (1,700 V).
##' - **Mass spectrometry**: Thermo Scientific Q-Exactive (MS1
##'   resolution = 70,000; MS1 accumulation time = 300ms; MS2
##'   resolution = 70,000).
##' - **Data analysis**: MaxQuant(2.4.13.0) + DART-ID.
##'
##' @section Data collection:
##'
##' The PSM data were collected from a shared Google Drive folder that
##' is accessible from the SlavovLab website (see `Source` section).
##' The folder ('/002-singleCellDataGeneration') contains the following
##' files of interest:
##'
##' - `ev_updated_NS.DIA.txt`: the MaxQuant/DART-ID output file
##' - `annotation.csv`: sample annotation
##' - `batch.csv`: batch annotation
##'
##' We combined the sample annotation and the batch annotation in
##' a single table. We also formatted the quantification table so that
##' columns match with those of the annotation and filter only for
##' single-cell runs. Both table are then combined in a single
##' [QFeatures] object using the [scp::readSCP()] function.
##'
##' The peptide data were taken from the same google drive folder
##' (`EpiToMesen.TGFB.nPoP_trial1_pepByCellMatrix_NSThreshDART_medIntCrNorm.txt`).
##' The data were formatted to a [SingleCellExperiment] object and the sample
##' metadata were matched to the column names (mapping is retrieved
##' after running the SCoPE2 R script, `EMTTGFB_singleCellProcessing.R`) and
##' stored in the `colData`. The object is then added to the [QFeatures] object
##' and the rows of the PSM data are linked to the rows of the peptide data
##' based on the peptide sequence information through an `AssayLink` object.
##'
##' The imputed protein data were taken from the same google drive folder
##' (`EpiToMesen.TGFB.nPoP_trial1_1PercDartFDRTMTBulkDIA.WallE_imputed.txt`).
##' The data were formatted to a [SingleCellExperiment] object and the sample
##' metadata were matched to the column names (mapping is retrieved
##' after running the SCoPE2 R script, `EMTTGFB_singleCellProcessing.R`) and
##' stored in the `colData`. The object is then added to the [QFeatures] object
##' and the rows of the peptide data are linked to the rows of the protein data
##' based on the protein sequence information through an `AssayLink` object.
##'
##' The unimputed protein data were taken from the same google drive folder
##' (`EpiToMesen.TGFB.nPoP_trial1_1PercDartFDRTMTBulkDIA.WallE_unimputed.txt`).
##' The data were formatted and added exactly as imputed data.
##'
##' @source
##' The data were downloaded from the
##' [Slavov Lab](https://scp.slavovlab.net/Khan_et_al_2023) website via a
##' shared Google Drive
##' [folder](https://drive.google.com/drive/folders/1zCsRKWNQuAz5msxx0DfjDrIe6pUjqQmj).
##' The raw data and the quantification data can also be found in the
##' MassIVE repository `MSV000092872`:
##' ftp://MSV000092872@massive.ucsd.edu/.
##'
##' @references
##' Saad Khan, Rachel Conover, Anand R. Asthagiri, Nikolai Slavov. 2023.
##' "Dynamics of single-cell protein covariation during epithelial–mesenchymal
##' transition." Journal of Proteome Research.
##' ([link to article](https://pubs.acs.org/doi/10.1021/acs.jproteome.4c00277)).
##'
##' @examples
##' \donttest{
##' khan2023()
##' }
##'
##' @keywords datasets
##'
"khan2023"
