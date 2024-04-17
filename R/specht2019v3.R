##' Specht et al. 2019 - SCoPE2 (biorRxiv): macrophages vs monocytes
##' (version 3)
##'
##' Single cell proteomics data acquired by the Slavov Lab. This is
##' the version 3 of the data released in October 2020. It contains
##' quantitative information of macrophages and monocytes at PSM,
##' peptide and protein level.
##'
##' @format A [QFeatures] object with 179 assays, each assay being a
##' [SingleCellExperiment] object:
##'
##' - Assay 1-63: PSM data for SCoPE2 sets acquired with a TMT-11plex
##'   protocol, hence those assays contain 11 columns. Columns
##'   hold quantitative information from single-cell channels, carrier
##'   channels, reference channels, empty (blank) channels and unused
##'   channels.
##' - Assay 64-177: PSM data for SCoPE2 sets acquired with a
##'   TMT-16plex protocol, hence those assays contain 16 columns.
##'   Columns hold quantitative information from single-cell channels,
##'   carrier channels, reference channels, empty (blank) channels and
##'   unused channels.
##' - `peptides`: peptide data containing quantitative data for 9208
##'   peptides and 1018 single-cells.
##' - `proteins`: protein data containing quantitative data for 2772
##'   proteins and 1018 single-cells.
##'
##' The `colData(specht2019v2())` contains cell type annotation and
##' batch annotation that are common to all assays. The description of
##' the `rowData` fields for the PSM data can be found in the
##' [`MaxQuant` documentation](http://www.coxdocs.org/doku.php?id=maxquant:table:evidencetable).
##'
##' @section Acquisition protocol:
##'
##' The data were acquired using the following setup. More information
##' can be found in the source article (see `References`).
##'
##' - **Cell isolation**: flow cytometry (BD FACSAria I).
##' - **Sample preparation** performed using the SCoPE2 protocol. mPOP
##'   cell lysis + trypsin digestion + TMT-11plex or 16plex labeling
##'   and pooling.
##' - **Separation**: online nLC (DionexUltiMate 3000 UHPLC with a
##'   25cm x 75um IonOpticksAurora Series UHPLC column; 200nL/min).
##' - **Ionization**: ESI (2,200V).
##' - **Mass spectrometry**: Thermo Scientific Q-Exactive (MS1
##'   resolution = 70,000; MS2 accumulation time = 300ms; MS2
##'   resolution = 70,000).
##' - **Data analysis**: DART-ID + MaxQuant (1.6.2.3).
##'
##' @section Data collection:
##'
##' The PSM data were collected from a shared Google Drive folder that
##' is accessible from the SlavovLab website (see `Source` section).
##' The folder contains the following files of interest:
##'
##' - `ev_updated_v2.txt`: the MaxQuant/DART-ID output file
##' - `annotation_fp60-97.csv`: sample annotation
##' - `batch_fp60-97.csv`: batch annotation
##'
##' We combined the sample annotation and the batch annotation in
##' a single table. We also formatted the quantification table so that
##' columns match with those of the annotation and filter only for
##' single-cell runs. Both table are then combined in a single
##' [QFeatures] object using the [scp::readSCP()] function.
##'
##' The peptide data were taken from the Slavov lab directly
##' (`Peptides-raw.csv`). It is provided as a spreadsheet. The data
##' were formatted to a [SingleCellExperiment] object and the sample
##' metadata were matched to the column names (mapping is retrieved
##' after running the SCoPE2 R script) and stored in the `colData`.
##' The object is then added to the [QFeatures] object (containing the
##' PSM assays) and the rows of the peptide data are linked to the
##' rows of the PSM data based on the peptide sequence information
##' through an `AssayLink` object.
##'
##' The protein data (`Proteins-processed.csv`) is formatted similarly
##' to the peptide data, and the rows of the proteins were mapped onto
##' the rows of the peptide data based on the protein sequence
##' information.
##'
##' @note Since version 2, a serious bug in the data were corrected
##' for TMT channels 12 to 16. Many more cells are therefore contained
##' in the data. Version 2 is maintained for backward compatibility.
##' Although the final version of the article was published in 2021,
##' we have kept `specht2019v3` as the data set name for consistency
##' with the previous data version `specht2019v2`.
##'
##' @source
##' The data were downloaded from the
##' [Slavov Lab](https://scope2.slavovlab.net/docs/data) website via a
##' shared Google Drive
##' [folder](https://drive.google.com/drive/folders/1VzBfmNxziRYqayx3SP-cOe2gu129Obgx).
##' The raw data and the quantification data can also be found in the
##' massIVE repository `MSV000083945`:
##' ftp://massive.ucsd.edu/MSV000083945.
##'
##' @references Specht, Harrison, Edward Emmott, Aleksandra A.
##' Petelski, R. Gray Huffman, David H. Perlman, Marco Serra, Peter
##' Kharchenko, Antonius Koller, and Nikolai Slavov. 2021.
##' "Single-Cell Proteomic and Transcriptomic Analysis of Macrophage
##' Heterogeneity Using SCoPE2." Genome Biology 22 (1): 50.
##' ([link to article](http://dx.doi.org/10.1186/s13059-021-02267-5)).
##'
##' @examples
##' \donttest{
##' specht2019v3()
##' }
##'
##' @keywords datasets
##'
"specht2019v3"
