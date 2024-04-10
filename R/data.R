####---- specht2019v2 ----####


##' Specht et al. 2019 - SCoPE2 (biorRxiv): macrophages vs monocytes
##' (version 2)
##'
##' @description
##'
##' Single cell proteomics data acquired by the Slavov Lab. This is
##' the version 2 of the data released in December 2019. It contains
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
##'   cell lysis + trypsin digestion + TMT-11plex or 16plex labelling
##'   and pooling.
##' - **Separation**: online nLC (DionexUltiMate 3000 UHPLC with a
##'   25cm x 75um IonOpticksAurora Series UHPLC column; 200nL/min).
##' - **Ionization**: ESI (2,200V).
##' - **Mass spectrometry**: Thermo Scientific Q-Exactive (MS1
##'   resolution = 70,000; MS1 accumulation time = 300ms; MS2
##'   resolution = 70,000).
##' - **Data analysis**: DART-ID + MaxQuant (1.6.2.3).
##'
##' @section Data collection:
##'
##' The PSM data were collected from a shared Google Drive folder that
##' is accessible from the SlavovLab website (see `Source` section).
##' The folder contains the following files of interest:
##'
##' - `ev_updated.txt`: the MaxQuant/DART-ID output file
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
##' @source
##' The data were downloaded from the
##' [Slavov Lab](https://scope2.slavovlab.net/mass-spec/data) website via a
##' shared Google Drive
##' [folder](https://drive.google.com/drive/folders/1VzBfmNxziRYqayx3SP-cOe2gu129Obgx).
##' The raw data and the quantification data can also be found in the
##' massIVE repository `MSV000083945`:
##' ftp://massive.ucsd.edu/MSV000083945.
##'
##' @references Specht, Harrison, Edward Emmott, Aleksandra A.
##' Petelski, R. Gray Huffman, David H. Perlman, Marco Serra, Peter
##' Kharchenko, Antonius Koller, and Nikolai Slavov. 2019.
##' "Single-Cell Mass-Spectrometry Quantifies the Emergence of
##' Macrophage Heterogeneity." bioRxiv.
##' ([link to article](https://doi.org/10.1101/665307)).
##'
##' @examples
##' \donttest{
##' specht2019v2()
##' }
##'
##' @keywords datasets
##'
"specht2019v2"

####---- specht2019v3 ----####

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


####---- dou2019_lysates ----####


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


####---- dou2019_mouse ----####


##' Dou et al. 2019 (Anal. Chem.): murine cell lines
##'
##' @description
##'
##' Single-cell proteomics using nanoPOTS combined with TMT isobaric
##' labeling. It contains quantitative information at PSM and protein
##' level. The cell types are either "Raw" (macrophage cells), "C10"
##' (epihelial cells), or "SVEC" (endothelial cells). Out of the 132
##' wells, 72 contain single cells, corresponding to 24 C10 cells, 24
##' RAW cells, and 24 SVEC. The other wells are either boosting
##' channels (12), empty channels (36) or reference channels (12).
##' Boosting and reference channels are balanced (1:1:1) mixes of C10,
##' SVEC, and RAW samples at 5 ng and 0.2 ng, respectively. The
##' different cell types where evenly distributed across 4 nanoPOTS
##' chips. Samples were 11-plexed with TMT labeling.
##'
##' @format A [QFeatures] object with 13 assays, each assay being a
##' [SingleCellExperiment] object:
##'
##' - `Single_Cell_Chip_X_Y`: PSM data with 11 columns corresponding
##'   to the TMT channels (see `Notes`). The `X` indicates the chip
##'   number (from 1 to 4) and `Y` indicates the row name on the chip
##'   (from A to C).
##' - `peptides`: peptide data containing quantitative data for 15,492
##'   peptides in 132 samples (run 1 and run 2 combined).
##' - `proteins`: protein data containing quantitative data for 2331
##'   proteins in 132 samples (all runs combined).
##'
##' Sample annotation is stored in `colData(dou2019_mouse())`.
##'
##' @section Acquisition protocol:
##'
##' The data were acquired using the following setup. More information
##' can be found in the source article (see `References`).
##'
##' - **Cell isolation**: single-cells from the three murine cell
##'   lines were isolated using FACS (BD Influx II cell sorter ).
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
##'
##' - `Single_Cell_Chip_*_*_msgfplus.mzid`: the MS-GF+ identification
##'   result files.
##' - `Single_Cell_Chip_*_*_ReporterIons.txt`: the MASIC
##'   quantification result files.
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
##' supplied as an Excel file `ac9b03349_si_005.xlsx`. The file
##' contains 7 sheets from which we only took the 2nd (named
##' `01 - Raw sc protein data`) with the combined protein data for the
##' 12 runs. We converted the data to a [SingleCellExperiment] object
##' and added the object as a new assay in the [QFeatures] dataset
##' (containing the PSM data). Links between the proteins and the
##' corresponding PSM were created. Note that links to protein groups
##' are not stored.
##'
##' @note Although a TMT-10plex labeling is reported in the article,
##' the PSM data contained 11 channels for each run. Those 11th
##' channel contain mostly missing data and are hence assumed to be
##' empty channels.
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
##' [dou2019_lysates], [dou2019_boosting]
##'
##' @examples
##' \donttest{
##' dou2019_mouse()
##' }
##'
##' @keywords datasets
##'
"dou2019_mouse"


####---- dou2019_boosting ----####


##' Dou et al. 2019 (Anal. Chem.): testing boosting ratios
##'
##' @description
##'
##' Single-cell proteomics using nanoPOTS combined with TMT isobaric
##' labeling. It contains quantitative information at PSM and protein
##' level. The cell types are either "Raw" (macrophage cells), "C10"
##' (epihelial cells), or "SVEC" (endothelial cells). Each cell is
##' replicated 2 or 3 times. Each cell type was run using 3 levels of
##' boosting: 0 ng (no boosting), 5 ng or 50 ng. When boosting was
##' applied, 1 reference well and 1 boosting well were added,
##' otherwise 1 empty well was added. Each boosting setting (0ng, 5ng,
##' 50ng) was run in duplicate.
##'
##' @format A [QFeatures] object with 7 assays, each assay being a
##' [SingleCellExperiment] object:
##'
##' - `Boosting_X_run_Y`: PSM data with 10 columns corresponding to
##'   the TMT-10plex channels. The `X` indicates the boosting amount
##'   (0ng, 5ng or 50ng) and `Y` indicates the run number (1 or 2).
##' - `peptides`: peptide data containing quantitative data for 13,462
##'   peptides in 60 samples (run 1 and run 2 combined).
##' - `proteins`: protein data containing quantitative data for 1436
##'   proteins and 60 samples (all runs combined).
##'
##' Sample annotation is stored in `colData(dou2019_boosting())`.
##'
##' @section Acquisition protocol:
##'
##' The data were acquired using the following setup. More information
##' can be found in the source article (see `References`).
##'
##' - **Cell isolation**: single-cells from the three murine cell
##'   lines were isolated using FACS (BD Influx II cell sorter ).
##'   Boosting sample were prepared (presumably in bulk) from 1:1:1
##'   mix of the three cell lines.
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
##' - `Boosting_*ng_run_*_msgfplus.mzid`: the MS-GF+ identification
##'   result files.
##' - `Boosting_*ng_run_*_ReporterIons.txt`: the MASIC quantification
##'   result files.
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
##' supplied as an Excel file `ac9b03349_si_004.xlsx`. The file
##' contains 7 sheets from which we took the 2nd, 4th and 6th sheets
##' (named `01 - No Boost raw data`, `03 - 5ng boost raw data`,
##' `05 - 50ng boost raw data`, respectively). The sheets contain the
##' combined protein data for the duplicate runs given the boosting
##' amount. We joined the data for all boosting ration based on the
##' protein name and converted the data to a [SingleCellExperiment]
##' object. We then added the object as a new assay in the [QFeatures]
##' dataset (containing the PSM data). Links between the proteins and
##' the corresponding PSM were created. Note that links to protein
##' groups are not stored.
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
##' [dou2019_lysates], [dou2019_mouse]
##'
##' @examples
##' \donttest{
##' dou2019_boosting()
##' }
##'
##' @keywords datasets
##'
##'
"dou2019_boosting"


####---- zhu2018MCP ----####


##' Zhu et al. 2018 (Mol. Cel. Prot.): rat brain laser dissections
##'
##' Near single-cell proteomics data of laser captured
##' micro-dissection samples. The samples are 24 brain sections from
##' rat pups (day 17). The slices are 12 um thick squares of either
##' 50, 100, or 200 um width. 5 samples were dissected from the corpus
##' callum (`CC`), 4 samples were dissected from the corpus collosum
##' (`CP`), 13 samples were extracted from the cerebral cortex
##' (`CTX`), and 2 samples are labeled as (`Mix`).
##'
##' @format A [QFeatures] object with 4 assays, each assay being a
##' [SingleCellExperiment] object:
##'
##' - `peptides`: quantitative information for 13,055 peptides from
##'   24 samples
##' - `proteins_intensity`: protein intensities for 2,257 proteins
##'   from 24 samples
##' - `proteins_LFQ`: LFQ intensities for 2,257 proteins from 24 samples
##' - `proteins_iBAQ`: iBAQ values for 2,257 proteins from 24 samples
##'
##' Sample annotation is stored in `colData(zhu2018MCP())`.
##'
##' @section Acquisition protocol:
##'
##' The data were acquired using the following setup. More information
##' can be found in the original article (see `References`).
##'
##' - **Cell isolation**: brain patches were collected using
##'   laser-capture microdissection (PALM MicroBeam) on flash frozen
##'   rat (*Rattus norvergicus*) brain tissues. Note that the samples
##'   were stained with H&E before dissection for histological
##'   analysis. DMSO is used as sample collection solution
##' - **Sample preparation** performed using the nanoPOTs device: DMSO
##'   evaporation + protein extraction (DMM + DTT) + alkylation (IAA)
##'   + Lys-C digestion + trypsin digestion.
##' - **Separation**: nanoLC (Dionex UltiMate with an in-house packed
##'   60cm x 30um LC columns; 50nL/min)
##' - **Ionization**: ESI (2,000V)
##' - **Mass spectrometry**: Thermo Fisher Orbitrap Fusion Lumos
##'   Tribrid (MS1 accumulation time = 246ms; MS1 resolution =
##'   120,000; MS1 AGC = 3E6). The MS/MS settings  depend on the
##'   sample size, excepted for the AGC = 1E5. 50um (time = 502ms;
##'   resolution = 240,000), 100um (time = 246ms; resolution =
##'   120,000), 200um (time = 118ms; resolution = 60,000).
##' - **Data analysis**: MaxQuant (v1.5.3.30) + Perseus (v1.5.6.0) +
##'   Origin Pro 2017
##'
##' @section Data collection:
##'
##' The data were collected from the PRIDE repository (accession
##' ID: PXD008844).  We downloaded the `MaxQuant_Peptides.txt`
##' and the `MaxQuant_ProteinGroups.txt` files containing the
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
##' PXD008844. FTP link
##' ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2018/07/PXD008844
##'
##' @references
##' Zhu, Ying, Maowei Dou, Paul D. Piehowski, Yiran Liang, Fangjun
##' Wang, Rosalie K. Chu, William B. Chrisler, et al. 2018. “Spatially
##' Resolved Proteome Mapping of Laser Capture Microdissected Tissue
##' with Automated Sample Transfer to Nanodroplets.” Molecular &
##' Cellular Proteomics: MCP 17 (9): 1864–74
##' ([link to article](http://dx.doi.org/10.1074/mcp.TIR118.000686)).
##'
##' @examples
##' \donttest{
##' zhu2018MCP()
##' }
##'
##' @keywords datasets
##'
##'
"zhu2018MCP"


####---- zhu2018NC_hela ----####


##' Zhu et al. 2018 (Nat. Comm.): HeLa titration
##'
##' Near single-cell proteomics data of HeLa samples containing
##' different number of cells. There are three groups of cell
##' concentrations: low (10-14 cells), medium (35-45 cells) and high
##' (137-141 cells). The data also contain measures for blanks, HeLa
##' lysates (50 cell equivalent) and 2 cancer cell line lysates (MCF7
##' and THP1, 50 cell equivalent).
##'
##' @format A [QFeatures] object with 4 assays, each assay being a
##' [SingleCellExperiment] object:
##'
##' - `peptides`: quantitative information for 37,795 peptides from
##'   21 samples
##' - `proteins_intensity`: protein intensities for 3,984 proteins
##'   from 21 samples
##' - `proteins_LFQ`: LFQ intensities for 3,984 proteins from 21
##'   samples
##' - `proteins_iBAQ`: iBAQ values for 3,984 proteins from 21 samples
##'
##' Sample annotation is stored in `colData(zhu2018NC_hela())`.
##'
##' @section Acquisition protocol:
##'
##' The data were acquired using the following setup. More information
##' can be found in the original article (see `References`).
##'
##' - **Cell isolation**: HeLa cell concentration was adjusted by
##'   serial dilution and cell counting was performed manually using
##'   an inverted microscope.
##' - **Sample preparation** performed using the nanoPOTs device.
##'   Protein extraction using RapiGest (+ DTT) + alkylation (IAA) +
##'   Lys-C digestion + cleave RapiGest (formic acid).
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
##'   2017
##'
##' @section Data collection:
##'
##' The data were collected from the PRIDE repository (accession
##' ID: PXD006847).  We downloaded the `CulturedCells_peptides.txt`
##' and the `CulturedCells_proteinGroups.txt` files containing the
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
##' PXD006847. FTP link:
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
##' @seealso The same experiment was conducted on HeLa lysates:
##' [zhu2018NC_lysates].
##'
##' @examples
##' \donttest{
##' zhu2018NC_hela()
##' }
##'
##' @keywords datasets
##'
"zhu2018NC_hela"


####---- zhu2018NC_lysates ----####


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


####---- zhu2018NC_islets ----####


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


####---- cong2020AC ----####


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


####---- zhu2019EL ----####


##' Zhu et al. 2019 (eLife): chicken utricle cells
##'
##'
##' Single-cell proteomics data from chicken utricle acquired to
##' study the hair-cell development. The cells are isolated from
##' peeled utrical epithelium and separated into hair cells (FM1-43
##' high) and supporting cells (FM1-43 low). The sample contain either
##' 1 cell (n = 28), 3 cells (n = 7), 5 cells (n = 8) or 20 cells (n =
##' 14).
##'
##' @format A [QFeatures] object with 62 assays, each assay being a
##' [SingleCellExperiment] object:
##'
##' - `XYZw`: 60 assays containing PSM data. The sample are annotated
##'   as follows. `X` indicates the experiment, either 1 or 2. `Y`
##'   indicated the FM1-43 signal, either high (H) or low (L). `Z`
##'   indicates the number of cells (0, 1, 3, 5 or 20). `w` indicates
##'   the replicate, starting from `a`, it can go up to `j`.
##' - `peptides`: quantitative data for 3444 peptides in 60 samples
##'   (all runs are combined).
##' - `proteins_intensity`: protein intensities for 840 proteins
##'   from 24 samples
##' - `proteins_iBAQ`: iBAQ values for 840 proteins from 24 samples
##'
##' Sample annotation is stored in `colData(zhu2019EL())`.
##'
##' @section Acquisition protocol:
##'
##' The data were acquired using the following setup. More information
##' can be found in the source article (see `References`).
##'
##' - **Cell isolation**: The cells were taken from the utricles of
##'   E15 chick embryos. Samples were stained with FM1-43FX and the
##'   cells were dissociated using enzymatic digestion. Cells were
##'   FACS sorted (BD Influx) and split based on their FM1-43 signal,
##'   while ensuring no debris, doublets or dead cells are retained.
##' - **Sample preparation** performed using the nanoPOTs device. Cell
##'   lysis and protein extraction and reduction are performed using
##'   dodecyl beta-D-maltoside + DTT + ammonium bicarbonate. Protein
##'   were then alkylated using IAA. Protein digestion is performed
##'   using Lys-C and trypsin. Finally samples acidification is
##'   performed using formic acid.
##' - **Separation**: Dionex UltiMate pump with an C18-Packed column
##'   (50cm x 30um; 60nL/min)
##' - **Ionization**: ESI (2,000V)
##' - **Mass spectrometry**: Orbitrap Fusion Lumos Tribrid. MS1
##'   settings: accumulation time = 246ms; resolution = 120,000; AGC =
##'   3E6. MS/MS settings: accumulation time = 502ms; resolution =
##'   120,000; AGC = 2E5.
##' - **Data analysis**: Andromeda & MaxQuant (v1.5.3.30) and the
##'   search database is NCBI GRCg6a.
##'
##' @section Data collection:
##'
##' All data were collected from the PRIDE repository (accession ID:
##' PXD014256).
##'
##' The sample annotation information is provided in the
##' `Zhu_2019_chick_single_cell_samples_CORRECTED.xlsx` file. This file
##' was given during a personal discussion and is a corrected version
##' of the annotation table available on the PRIDE repository.
##'
##' The PSM data were found in the `evidence.txt` (in the
##' `Experiment 1+ 2`) folder. The PSM data were filtered so that it
##' contains only samples that are annotated. The data were then
##' converted to a [QFeatures] object using the [scp::readSCP()]
##' function.
##'
##' The peptide data were found in the `peptides.txt` file. The column
##' names holding the quantitative data were adapted to match the
##' sample names in the [QFeatures] object. The data were then
##' converted to a [SingleCellExperiment] object and then inserted in
##' the [QFeatures] object. Links between the PSMs and the peptides
##' were added
##'
##' A similar procedure was applied to the protein data. The data were
##' found in the `proteinGroups.txt` file. We split the protein table
##' to separate the two types of quantification: summed intensity and
##' intensity based absolute quantification (iBAQ). Both tables are
##' converted to [SingleCellExperiment] objects and are added to the
##' [QFeatures] object as well as the `AssayLink` between peptides and
##' proteins.
##'
##' @source
##' The PSM data can be downloaded from the PRIDE repository
##' PXD014256. The source link is:
##' ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/11/PXD014256
##'
##' @references
##'
##' Zhu, Ying, Mirko Scheibinger, Daniel Christian Ellwanger, Jocelyn
##' F. Krey, Dongseok Choi, Ryan T. Kelly, Stefan Heller, and Peter G.
##' Barr-Gillespie. 2019. “Single-Cell Proteomics Reveals Changes in
##' Expression during Hair-Cell Development.” eLife 8 (November).
##' ([link to article](https://doi.org/10.7554/eLife.50777)).
##'
##' @examples
##' \donttest{
##' zhu2019EL()
##' }
##'
##' @keywords datasets
##'
"zhu2019EL"


####---- liang2020_hela ----####


##' Liang et al. 2020 (Anal. Chem.): HeLa cells (MaxQuant preprocessing)
##'
##' Single-cell proteomics data from HeLa cells using the autoPOTS
##' acquisition workflow. The samples contain either no cells (blanks),
##' 1 cell, 10 cells, 150 cells or 500 cells. Samples containing
##' between 0 and 10 cells are isolated using micro-pipetting while
##' samples containing between 150 and 500 cells were prepared using
##' dilution of a bulk sample.
##'
##' @format A [QFeatures] object with 17 assays, each assay being a
##' [SingleCellExperiment] object:
##'
##' - `HeLa_*`: 15 assays containing PSM data.
##' - `peptides`: quantitative data for 48705 peptides in 15 samples
##'   (all runs are combined).
##' - `proteins`: quantitative data for 3970 protein groups in 15
##'   samples (all runs combined).
##'
##' Sample annotation is stored in `colData(liang2020_hela())`.
##'
##' @section Acquisition protocol:
##'
##' The data were acquired using the following setup. More information
##' can be found in the source article (see `References`).
##'
##' - **Cell isolation**: The HeLa cells come from a commercially
##'   available cell line. Samples containing between 0 and 10 cells
##'   were isolated using micro-manipulation and the counts were
##'   validated using a microscope. Samples containing between 150 and
##'   500 cells were prepared by diluting a bulk sample and the exact
##'   counts were evaluated by obtaining phtotmicrographs.
##' - **Sample preparation** performed using the autoPOTS worflow that
##'   relied on the OT-2 pipeting robot. Cell are lysed using
##'   sonication. Samples are then processed by successive incubation
##'   with DTT (reduction), then IAA (alkylation), then Lys-C and
##'   trypsin (protein digestion).
##' - **Separation**: Samples were injected on the column using a
##'   modified Ultimate WPS-3000 TPL autosampler coupled to an UltiMate
##'   3000 RSLCnano pump. The LC column is a home-packed nanoLC column
##'   (45cm x 30um; 40nL/min)
##' - **Ionization**: Nanospray Flex ion source (2,000V)
##' - **Mass spectrometry**: Orbitrap Exploris 480. MS1 settings:
##'   accumulation time = 250 ms (0-10 cells) or 100 ms (150-500 cells);
##'   resolution = 120,000; AGC = 100\%. MS2 settings: exlusion
##'   duration = 90 s (0-10 cells) or 60 s (150-500 cells) ; accumulation
##'   time = 500 ms (0-1 cell), 250 ms (10 cells), 100 ms (150 cells)
##'   or 50 ms (500 cells); resolution = 60,000 (0-10 cells) or 30,000
##'   (150-500 cells); AGC = 5E3 (0-1 cells) or 1E4 (10-500 cells).
##' - **Data analysis**: MaxQuant (v1.6.7.0) and the search database
##'   is Swiss-Prot (July 2020).
##'
##' @section Data collection:
##'
##' All data were collected from the PRIDE repository (accession ID:
##' PXD021882).
##'
##' The sample annotations were collected from the methods section and
##' from table S3 in the paper.
##'
##' The PSM data were found in the `evidence.txt` file. The data were
##' converted to a [QFeatures] object using the [scp::readSCP()]
##' function.
##'
##' The peptide data were found in the `peptides.txt` file. The column
##' names holding the quantitative data were adapted to match the
##' sample names in the [QFeatures] object. The data were then
##' converted to a [SingleCellExperiment] object and then inserted in
##' the [QFeatures] object. Links between the PSMs and the peptides
##' were added
##'
##' A similar procedure was applied to the protein data. The data were
##' found in the `proteinGroups.txt` file. The column names were
##' adapted, the data were converted to a [SingleCellExperiment]
##' object and then inserted in the [QFeatures] object. Links between
##' the peptides and the proteins were added
##'
##' @source
##' The PSM data can be downloaded from the PRIDE repository
##' PXD021882 The source link is:
##' http://ftp.pride.ebi.ac.uk/pride/data/archive/2020/12/PXD021882/
##'
##' @references
##'
##' Liang, Yiran, Hayden Acor, Michaela A. McCown, Andikan J. Nwosu,
##' Hannah Boekweg, Nathaniel B. Axtell, Thy Truong, Yongzheng Cong,
##' Samuel H. Payne, and Ryan T. Kelly. 2020. “Fully Automated Sample
##' Processing and Analysis Workflow for Low-Input Proteome
##' Profiling.” Analytical Chemistry, December.
##' ([link to article](https://doi.org/10.1021/acs.analchem.0c04240)).
##'
##' @examples
##' \donttest{
##' liang2020_hela()
##' }
##'
##' @keywords datasets
##'
"liang2020_hela"



####---- schoof2021 ----####


##' Schoof et al. 2021 (Nat. Comm.): acute myeloid leukemia
##' differentiation
##'
##' Single-cell proteomics data from OCI-AML8227 cell culture to
##' reconstruct the cellular hierarchy. The data were acquired using
##' TMTpro multiplexing. The samples contain either no cells,
##' single cells, 10 cells (reference channel) 200 cells (booster
##' channel) or are simply empty wells. Single cells are expected to
##' be one of progenitor cells (`PROG`), leukaemia stem cells (`LSC`),
##' CD38- blast cells (`BLAST CD38-`) or CD38+ blast cells
##' (`BLAST CD38+`). Booster are either a known 1:1:1 mix of cells
##' (PROG, LSC and BLAST) or are isolated directly from the bulk
##' sample. Samples were isolated and annotated using flow cytometry.
##'
##' @format A [QFeatures] object with 194 assays, each assay being a
##' [SingleCellExperiment] object:
##'
##' - `F*`: 192 assays containing PSM quantification data for 16
##'    TMT channels. The quantification data contain signal to noise
##'    ratios as computed by Proteome Discoverer.
##' - `proteins`: quantitative data for 2898 protein groups in 3072
##'   samples (all runs combined). The quantification data contain
##'   signal to noise ratios as computed by Proteome Discoverer.
##' - `logNormProteins`: quantitative data for 2723 protein groups in
##'   2025 single-cell samples. This assay is the protein datasets that
##'   was processed by the authors. Dimension reduction and clustering
##'   data are also available in the `reducedDims` and `colData` slots,
##'   respectively
##'
##' Sample annotation is stored in `colData(schoof2021())`. The cell
##' type annotation is stored in the `Population` column. The flow
##' cytometry data is also available: FSC-A, FSC-H, FSC-W, SSC-A,
##' SSC-H, SSC-W, APC-Cy7-A (= CD34) and PE-A (= CD38).
##'
##' @section Acquisition protocol:
##'
##' The data were acquired using the following setup. More information
##' can be found in the source article (see `References`).
##'
##' - **Sample isolation**: cultured AML 8227 cells were stained with
##'   anti-CD34 and anti-CD38. The sorting was performed by FACSAria
##'   instrument and deposited in 384 well plates.
##' - **Sample preparation**: cells are lysed using freeze-boil and
##'   sonication in a lysis buffer (TFE) that also includes reduction
##'   and alkylation reagents (TCEP and CAA), followed by trypsin
##'   (protein) and benzonase (DNA) digestion, TMT-16 labeling and
##'   quenching, desalting using SOLAµ C18 plate, peptide
##'   concentration, pooling and peptide concentration again. The
##'   booster channel contains 200 cell equivalents.
##' - **Liquid chromatography**: peptides are separated using a C18
##'   reverse-phase column (50cm x 75 µm i.d., Thermo EasySpray) combined
##'   to a Thermo EasyLC 1200 for 160 minute gradient with a flowrate of
##'   100nl/min.
##' - **Mass spectrometry**: FAIMSPro interface is used. MS1 setup:
##'   resolution 60.000, AGC target of 300%, accumulation of 50ms. MS2
##'   setup: resolution 45.000, AGC target of 150, 300 or 500%,
##'   accumulation of 150, 300, 500, or 1000ms.
##' - **Raw data processing**: Proteome Discoverer 2.4 + Sequest spectral
##'   search engine and validation with Percolator
##'
##' @section Data collection:
##'
##' All data were collected from the PRIDE repository (accession ID:
##' PXD020586). The data and metadata were extracted from the
##' `SCeptre_FINAL.zip` file.
##'
##' We performed extensive data wrangling to combine al the metadata
##' available from different files into a single table available using
##' `colData(schoof2021)`.
##'
##' The PSM data were found in the `bulk_PSMs.txt` file. Contaminants
##' were defined based on the protein accessions listed in
##' `contaminant.txt`. The data were converted to a [QFeatures]
##'  object using the [scp::readSCP()] function.
##'
##' The protein data were found in the `bulk_Proteins.txt` file.
##' Contaminants were defined based on the protein accessions listed
##' in `contaminant.txt`.The column names holding the quantitative
##' data were adapted to match the sample names in the [QFeatures]
##' object. Unnecessary feature annotations (such as in which assay
##' a protein is found) were removed. Feature names were created
##' following the procedure in SCeptre: features names are the
##' protein symbol (or accession if missing) and if duplicated
##' symbols are present (protein isoforms), they are made unique by
##' appending the protein accession.  Contaminants were defined based
##' on the protein accessions listed in `contaminant.txt`. The data
##' were then converted to a [SingleCellExperiment] object and
##' inserted in the [QFeatures] object.
##'
##' The log-normalized protein data were found in the `bulk.h5ad` file.
##' This dataset was generated by the authors by running the notebook
##' called `bulk.ipynb`. The `bulk.h5ad` was loaded as an `AnnData`
##' object using the `scanpy` Python module. The object was then
##' converted to a `SingleCellExperiment` object using the
##' `zellkonverter` package. The column names holding the quantitative
##' data were adapted to match the sample names in the [QFeatures]
##' object. The data were then inserted in the [QFeatures] object.
##'
##' The script to reproduce the `QFeatures` object is available at
##' `system.file("scripts", "make-data_schoof2021.R", package = "scpdata")`
##'
##' @source
##'
##' The PSM and protein data can be downloaded from the PRIDE
##' repository PXD020586 The source link is:
##' https://www.ebi.ac.uk/pride/archive/projects/PXD020586
##'
##' @references
##'
##' Schoof, Erwin M., Benjamin Furtwängler, Nil Üresin, Nicolas Rapin,
##' Simonas Savickas, Coline Gentil, Eric Lechman, Ulrich auf Dem
##' Keller, John E. Dick, and Bo T. Porse. 2021. “Quantitative
##' Single-Cell Proteomics as a Tool to Characterize Cellular
##' Hierarchies.” Nature Communications 12 (1): 745679.
##' ([link to article](http://dx.doi.org/10.1038/s41467-021-23667-y)).
##'
##' @examples
##' \donttest{
##' schoof2021()
##' }
##'
##' @keywords datasets
##'
"schoof2021"


####---- williams2020 LFQ ----####


##' Williams et al. 2020 (Anal. Chem.): MCF10A cell line
##'
##' Single-cell label free proteomics data from a MCF10A cell line
##' culture. The data were acquired using a label-free quantification
##' protocole based on the nanoPOTS technology. The objective was to
##' test 2 elution gradients for single-cell applications and to
##' demonstrate successful use of the new nanoPOTS autosampler
##' presented in the article. The samples contain either no cells,
##' single cells, 3 cells, 10 cells  50 cells.
##'
##' @format A [QFeatures] object with 9 assays, each assay being a
##' [SingleCellExperiment] object:
##'
##' - `peptides_[30 or 60]min_[intensity or LFQ]`: 3 assays
##'   containing peptide intensities or LFQ normalized
##'   quantifications (see `References`) for either a 30min or a 60 min
##'   gradient.
##' - `proteins_[30 or 60]min_[intensity or iBAQ or LFQ]`: 6 assays
##'   containing protein intensities, iBAQ normalized or LFQ normalized
##'   quantifications (see `References`) for either a 30min or a 60 min
##'   gradient.
##'
##' Sample annotation is stored in `colData(williams2020_lfq())`.
##'
##' @section Acquisition protocol:
##'
##' The data were acquired using the following setup. More information
##' can be found in the source article (see `References`).
##'
##' - **Sample isolation**: cultured MCF10A cells were isolated using
##'   flow-cytometry based cell sorting and deposit on nanoPOTS
##'   microwells
##' - **Sample preparation**: cells are lysed using using a DDM+DTT
##'   lysis buffer. Alkylation was then performed using an IAA solution.
##'   Proteins are digested with Lys-C and trypsin followed by
##'   acidification with FA. Sample droplets are then dried until
##'   LC-MS/MS analysis.
##' - **Liquid chromatography**: peptides are loaded using the new
##'   autosampler described in the paper. Samples are loaded using a
##'   a homemade miniature syringe pump. The samples are then desalted
##'   and concentrated through a SPE column (4cm x 100µm i.d. packed
##'   with 5µm C18) with microflow LC pump. The peptides are then eluted
##'   from a long LC column (60cm x 50 µm i.d. packed with 3µm C18)
##'   coupled to a nanoflox LC pump at 150nL/mL with either a 30 min
##'   or a 60 min gradient.
##' - **Mass spectrometry**: MS/MS was performed on an Orbitrap Fusion
##'   Lumos Tribrid MS coupled to a 2kV ESI. MS1 setup: Orbitrap
##'   analyzer at resolution 120.000, AGC target of 1E6, accumulation
##'   of 246ms. MS2 setup: ion trap with CID at resolution 60.000, AGC
##'   target of 2E4, accumulation of 120ms (50 cells) or 250ms (0-10
##'   cells).
##' - **Raw data processing**: preprocessing using Maxquant v1.6.2.10
##'   that use Andromeda search engine (with UniProtKB 2016-21-29),
##'   MBR and LFQ normalization were enabled.
##'
##' @section Data collection:
##'
##' All data were collected from the MASSIVE repository (accession ID:
##' MSV000085230).
##'
##' The peptide and protein data were extracted from the `Peptides_[...].txt`
##' or `ProteinGroups[...].txt` files, respectively, in the
##' `MCF10A_LC_[30 or 60]minutes` folders.
##'
##' The tables were duplicated so that peptide intensisities, peptide
##' LFQ, protein intensities, protein LFQ and protein intensities are
##' contained in separate tables. Tables are then converted to
##' [SingleCellExperiment] objects. Sample annotations were infered
##' from the sample names and from the paper. All data is combined in
##' a [QFeatures] object. [AssayLinks] were stored between peptide
##' assays and their corresponding proteins assays based on the
##' leading razor protein (hence only unique peptides are linked to
##' proteins).
##'
##' The script to reproduce the `QFeatures` object is available at
##' `system.file("scripts", "make-data_williams2020_lfq.R", package = "scpdata")`
##'
##' @section Suggestion:
##'
##' See `QFeatures::joinAssays` if you want to join the 30min and
##' 60min assays in a single assay for an integrated analysis.
##'
##' @source
##'
##' The PSM and protein data can be downloaded from the MASSIVE
##' repository MSV000085230.
##'
##' @references
##'
##' **Source article**: Williams, Sarah M., Andrey V. Liyu, Chia-Feng
##' Tsai, Ronald J. Moore, Daniel J. Orton, William B. Chrisler,
##' Matthew J. Gaffrey, et al. 2020. “Automated Coupling of
##' Nanodroplet Sample Preparation with Liquid Chromatography-Mass
##' Spectrometry for High-Throughput Single-Cell Proteomics.”
##' Analytical Chemistry 92 (15): 10588–96.
##' ([link to article](http://dx.doi.org/10.1021/acs.analchem.0c01551)).
##'
##' **LFQ normalization**: Cox, Jürgen, Marco Y. Hein, Christian A. Luber,
##' Igor Paron, Nagarjuna Nagaraj, and Matthias Mann. 2014. “Accurate
##' Proteome-Wide Label-Free Quantification by Delayed Normalization
##' and Maximal Peptide Ratio Extraction, Termed MaxLFQ.” Molecular
##' & Cellular Proteomics: MCP 13 (9): 2513–26.
##' ([link to article](http://dx.doi.org/10.1074/mcp.M113.031591)).
##'
##' **iBAQ normalization**: Schwanhäusser, Björn, Dorothea Busse, Na
##' Li, Gunnar Dittmar, Johannes Schuchhardt, Jana Wolf, Wei Chen, and
##' Matthias Selbach. 2011. “Global Quantification of Mammalian Gene
##' Expression Control.” Nature 473 (7347): 337–42.
##' ([link to article](http://dx.doi.org/10.1038/nature10098)).
##'
##' @examples
##' \donttest{
##' williams2020_lfq()
##' }
##'
##' @keywords datasets
##'
"williams2020_lfq"

####---- williams2020 TMT ----####


##' Williams et al. 2020 (Anal. Chem.): 3 AML cell line
##'
##' Single-cell label data from three acute myeloid
##' leukemia cell line culture (MOLM-14, K562, CMK). The data were
##' acquired using a TMT-based quantification protocole and the
##' nanoPOTS technology. The objective was to demonstrate successful
##' use of the new nanoPOTS autosampler presented in the source
##' article. The samples contain either carrier (10 ng), reference
##' (0.2ng), empty or single-cell samples..
##'
##' @format A [QFeatures] object with 4 assays, each assay being a
##' [SingleCellExperiment] object:
##'
##' - `peptides_[intensity or corrected]`: 2 assays containing peptide
##'   reporter ion intensities or corrected reporter ion intensities
##'   as computed by MaxQuant.
##' - `proteins_[intensity or corrected]`: 2 assays containing protein
##'   reporter ion intensities or corrected reporter ion intensities
##'   as computed by MaxQuant.
##'
##' Sample annotation is stored in `colData(williams2020_tmt())`.
##'
##' @section Acquisition protocol:
##'
##' The data were acquired using the following setup. More information
##' can be found in the source article (see `References`).
##'
##' - **Sample isolation**: cultured MOLM-14, K562 or CMK cells were
##'   isolated using flow-cytometry based cell sorting and deposit on
##'   nanoPOTS microwells
##' - **Sample preparation**: cells are lysed using using a DDM lysis
##'   buffer. Proteins are digested with trypsin followed by TMT
##'   labelling and quanching with HA. The samples are then acidified
##'   with FA, pooled in a single samples (adding carrier and reference
##'   peptide mixtures), and dried until LC-MS/MS analysis.
##' - **Liquid chromatography**: peptides are loaded using the new
##'   autosampler described in the paper. Samples are loaded using a
##'   a homemade miniature syringe pump. The samples are then desalted
##'   and concentrated through a SPE column (4cm x 100µm i.d. packed
##'   with 5µm C18) with microflow LC pump. The peptides are then eluted
##'   from a long LC column (60cm x 50 µm i.d. packed with 3µm C18)
##'   coupled to a nanoflox LC pump at 150nL/mL (elution time is not
##'   expliceted).
##' - **Mass spectrometry**: MS/MS was performed on an Orbitrap Fusion
##'   Lumos Tribrid MS coupled to a 2kV ESI. MS1 setup: Orbitrap
##'   analyzer at resolution 120.000, AGC target of 1E6, accumulation
##'   of 246ms. MS2 setup: Orbitrap with HCD at resolution 120.000, AGC
##'   target of 1E6, accumulation of 246ms.
##' - **Raw data processing**: preprocessing using Maxquant v1.6.2.10
##'   that use Andromeda search engine (with UniProtKB 2016-21-29).
##'
##' @section Data collection:
##'
##' All data were collected from the MASSIVE repository (accession ID:
##' MSV000085230).
##'
##' The peptide and protein data were extracted from the
##' `Peptides_AML_SingleCell.txt` or `ProteinGroups_AML_SingleCell.txt`
##' files, respectively, in the `AML_SingleCell` folders.
##'
##' The tables were duplicated so that intensisities and corrected
##' intensities are contained in separate tables. Tables are then
##' converted to [SingleCellExperiment] objects. Sample annotations
##' were inferred from the sample names, from table S2 and from the
##' Experimental Section of the paper. All data is combined in
##' a [QFeatures] object. [AssayLinks] were stored between peptide
##' assays and their corresponding proteins assays based on the
##' leading razor protein (hence only unique peptides are linked to
##' proteins).
##'
##' The script to reproduce the `QFeatures` object is available at
##' `system.file("scripts", "make-data_williams2020_tmt.R", package = "scpdata")`
##'
##' @source
##'
##' The PSM and protein data can be downloaded from the MASSIVE
##' repository MSV000085230.
##'
##' @references
##'
##' **Source article**: Williams, Sarah M., Andrey V. Liyu, Chia-Feng
##' Tsai, Ronald J. Moore, Daniel J. Orton, William B. Chrisler,
##' Matthew J. Gaffrey, et al. 2020. “Automated Coupling of
##' Nanodroplet Sample Preparation with Liquid Chromatography-Mass
##' Spectrometry for High-Throughput Single-Cell Proteomics.”
##' Analytical Chemistry 92 (15): 10588–96.
##' ([link to article](http://dx.doi.org/10.1021/acs.analchem.0c01551)).
##'
##' @examples
##' \donttest{
##' williams2020_tmt()
##' }
##'
##' @keywords datasets
##'
"williams2020_tmt"

####---- leduc2022_pSCoPE ----####

##' Leduc et al. 2022 - pSCoPE (biorRxiv): melanoma cells vs monocytes
##'
##' Single cell proteomics data acquired by the Slavov Lab. This is
##' the dataset associated to the third version of the preprint. It
##' contains quantitative information of melanoma cells and monocytes
##' at PSM, peptide and protein level. This version of the data was
##' acquired using the pSCoPE MS acquisition approach.
##'
##' @format A [QFeatures] object with 138 assays, each assay being a
##' [SingleCellExperiment] object:
##'
##' - Assay 1-134: PSM data acquired with a TMT-18plex protocol, hence
##'   those assays contain 18 columns. Columns hold quantitative
##'   information from single-cell channels, carrier channels,
##'   reference channels, empty (negative control) channels and
##'   unused channels.
##' - `peptides`: peptide data containing quantitative data for 20,804
##'   peptides and 1556 single-cells. These data have been filtered
##'   to keep high-quality PSMs, all batches have been normalized to
##'   the reference channel, PSMs were aggregated to peptides, and
##'   single-cells with low median coefficient of variation were kept.
##' - `peptides_log`: peptide data containing quantitative data for
##'   12,284 peptides and 1543 single-cells. The `peptides` data was
##'   further normalized, highly missing peptides were removed and the
##'   quantifications were log-transformed.
##' - `proteins_norm2`: protein data containing quantitative data for
##'   2844 proteins and 1543 single-cells. The peptides from
##'   `peptides_log` were aggregated to proteins and normalized.
##' - `proteins_processed`: protein data containing quantitative data
##'   for 2844 proteins and 1543 single-cells. The `proteins_norm2`
##'   data were imputed, batch corrected and normalized.
##'
##' The `colData(leduc2022_pSCoPE())` contains cell type annotation,
##' LC batch information, the TMT label, the MS run ID. We also added
##' the sample prep annotations provided by the cellenONE dispensing
##' device (only for single cells): time stamp of cell isolation by the
##' device, the diameter and elongation of the cell, the ID of the
##' sample glass side (4 slides in total), the field within the glass
##' (each slide is divided in 4 field), the pooled well ID (each field
##' contains 9 pools), the x and y coordinates of each cell dropped in
##' a field and of each cell pool upon pickup. Finally, we also
##' retrieved the melanoma subpopulation generated by the authors upon
##' data analysis. The main population is encoded as `A` while the
##' small population is encoded `B`. The description of the `rowData`
##' fields for the PSM data can be found in the
##' [`MaxQuant` documentation](http://www.coxdocs.org/doku.php?id=maxquant:table:evidencetable).
##'
##' @section Acquisition protocol:
##'
##' The data were acquired using the following setup. More information
##' can be found in the source article (see `References`).
##'
##' - **Cell isolation**: CellenONE cell sorting.
##' - **Sample preparation** performed using the improved SCoPE2
##'   protocol using the CellenONE liquid handling system. nPOP cell
##'   lysis (DMSO) + trypsin digestion + TMT-18plex
##'   labeling and pooling. A target library was generated as well to
##'   perform prioritized DDA (Huffman et al. 2022) using MaxQuant.Live
##'   (2.0.3).
##' - **Separation**: online nLC (DionexUltiMate 3000 UHPLC with a
##'   25cm x 75um IonOpticks Aurora Series UHPLC column; 200nL/min).
##' - **Ionization**: ESI (1,800V).
##' - **Mass spectrometry**: Thermo Scientific Q-Exactive (MS1
##'   resolution = 70,000; MS2 accumulation time = 300ms; MS2
##'   resolution = 70,000). Prioritized data acquisition was performed
##'   using the pSCoPE protocol (Huffman et al. 2022)
##' - **Data analysis**: MaxQuant (1.6.17.0) + DART-ID
##'
##' @section Data collection:
##'
##' The PSM data were collected from a shared Google Drive folder that
##' is accessible from the SlavovLab website (see `Source` section).
##' The folder contains the following files of interest:
##'
##' - `ev_updated.txt`: the MaxQuant/DART-ID output file
##' - `annotation.csv`: sample annotation
##' - `batch.csv`: batch annotation
##' - `t0.csv`: the processed data table containing the `peptides` data
##' - `t3.csv`: the processed data table containing the `peptides_log`
##'   data
##' - `t4b.csv`: the processed data table containing the
##'   `proteins_norm2` data
##' - `t6.csv`: the processed data table containing the
##'   `proteins_processed` data
##'
##' We combined the sample annotation and the batch annotation in
##' a single table. We also formatted the quantification table so that
##' columns match with those of the annotations. Both annotation and
##' quantification tables are then combined in a single [QFeatures]
##' object using the [scp::readSCP()] function.
##'
##' The 4 CSV files were loaded and formatted as [SingleCellExperiment]
##' objects and the sample metadata were matched to the column names
##' (mapping is retrieved after running the author's original R script)
##' and stored in the `colData`.
##' The object is then added to the [QFeatures] object (containing the
##' PSM assays) and the rows of the peptide data are linked to the
##' rows of the PSM data based on the peptide sequence information
##' through an `AssayLink` object.
##'
##' @source
##' The data were downloaded from the
##' [Slavov Lab](https://scp.slavovlab.net/Leduc_et_al_2022) website.
##' The raw data and the quantification data can also be found in the
##' massIVE repository `MSV000089159`:
##' ftp://massive.ucsd.edu/MSV000089159.
##'
##' @references
##' Andrew Leduc, Gray Huffman, and Nikolai Slavov. 2022. “Droplet
##' Sample Preparation for Single-Cell Proteomics Applied to the Cell
##' Cycle.” bioRxiv. [Link to article](https://doi.org/10.1101/2021.04.24.441211)
##'
##' Gray Huffman, Andrew Leduc, Christoph Wichmann, Marco di Gioia,
##' Francesco Borriello, Harrison Specht, Jason Derks, et al. 2022.
##' “Prioritized Single-Cell Proteomics Reveals Molecular and
##' Functional Polarization across Primary Macrophages.” bioRxiv.
##' [Link to article](https://doi.org/10.1101/2022.03.16.484655).
##'
##' @seealso
##' [leduc2022_plexDIA]
##'
##' @examples
##' \donttest{
##' leduc2022_pSCoPE()
##' }
##'
##' @keywords datasets
##'
"leduc2022_pSCoPE"

####---- leduc2022_plexDIA ----####

##' Leduc et al. 2022 - plexDIA (biorRxiv): melanoma cells
##'
##' Single cell proteomics data acquired by the Slavov Lab. This is
##' the dataset associated to the fourth version of the preprint (and
##' the Genome Biology publication). It contains quantitative
##' information of melanoma cells at precursor, peptide and protein level.
##' This version of the data was acquired using the plexDIA MS
##' acquisition protocol.
##'
##' @format A [QFeatures] object with 48 assays, each assay being a
##' [SingleCellExperiment] object:
##'
##' - Assay 1-45: precursor data acquired with a mTRAQ-3 protocol,
##'   hence those assays contain 3 columns. Columns hold quantitative
##'   information from single cells or negative control samples.
##' - `Ms1Extracted`: the DIA-NN MS1 extracted signal, it combines the
##'   information from assays 1-45.
##' - `peptides`: peptide data containing quantitative data for 3,608
##'   peptides and 104 single cells. The  data were filtered  to  1%
##'   protein FDR.
##' - `proteins`: protein data containing quantitative data for 508
##'   proteins and 105 single cells. Note that the peptide and protein
##'   data provided by the authors differ by 3 samples. The precursor
##'   data were aggregated to protein intensity using maxLFQ. The
##'   protein data were further median normalized by column and by row,
##'   log2 transformed, impute using KNN (k = 3), again median
##'   normalized by column and by row, batch corrected using ComBat,
##'   and median normalized by column and by row once more.
##'
##' The `colData(leduc2022_plexDIA())` contains cell type annotation and
##' batch annotation that are common to all assays. The description of
##' the `rowData` fields for the precursor data can be found in the
##' [`DIA-NN` documentation](https://github.com/vdemichev/DiaNN#readme).
##'
##' @section Acquisition protocol:
##'
##' The data were acquired using the following setup. More information
##' can be found in the source article (see `References`).
##'
##' - **Cell isolation**: CellenONE cell sorting.
##' - **Sample preparation** performed using the improved SCoPE2
##'   protocol using the CellenONE liquid handling system. nPOP cell
##'   lysis (DMSO) + trypsin digestion + mTRAQ-3
##'   labeling and pooling.
##' - **Separation**: online nLC (DionexUltiMate 3000 UHPLC with a
##'   25cm x 75um IonOpticks Aurora Series UHPLC column; 200nL/min).
##' - **Ionization**: ESI (1,800V).
##' - **Mass spectrometry**: Thermo Scientific Q-Exactive. The duty
##'   cycle = 1 MS1 + 4 DIA MS2 windows (120 Th, 120 Th, 200 Th and
##'   580 Th, spanning 378-1,402 m/z). Each MS1 and MS2 scan was
##'   conducted at 70,000 resolving power, 3×10E6 AGC and 300ms
##'   maximum injection time.
##' - **Data analysis**: DIA-NN.
##'
##' @section Data collection:
##'
##' The PSM data were collected from a shared Google Drive folder that
##' is accessible from the SlavovLab website (see `Source` section).
##' The folder contains the following files of interest:
##'
##' - `annotation_plexDIA.csv`: sample annotation
##' - `report_plexDIA_mel_nPOP.tsv`: the DIA-NN output file
##'   with the precursor data
##' - `report.pr_matrix_channels_ms1_extracted.tsv`: the DIA-NN
##'   output file with the combined precursor data
##' - `plexDIA_peptide.csv`: the processed data table containing the
##'   `peptide` data
##' - `plexDIA_protein_imputed.csv`: the processed data table
##'   containing the `protein` data
##'
##' We removed the failed runs as identified by the authors. We also
##' formatted the annotation and precuror quantification tables to
##' facilitate matching between corresponding columns. Both annotation
##' and quantification tables are then combined in a single [QFeatures]
##' object using `scp::readSCPfromDIANN()`.
##'
##' The `plexDIA_peptide.csv` and `plexDIA_protein_imputed.csv` files
##' were loaded and formatted as [SingleCellExperiment] objects. The
##' columns names were adapted to match those in the `QFeatures`
##' object. The `SingleCellExperiment` objects were then added to the
##' [QFeatures] object and the rows of the peptide data are linked to
##' the rows of the precursor data based on the peptide sequence or
##' the protein name through an `AssayLink` object.
##'
##' @source
##' The links to the data were found on the
##' [Slavov Lab website](https://scp.slavovlab.net/Leduc_et_al_2022).
##' The data were downloaded from the
##' [Google drive folder 1](https://drive.google.com/drive/folders/117ZUG5aFIJt0vrqIxpKXQJorNtekO-BV) and
##' [Google drive folder 2](https://drive.google.com/drive/folders/12-H2a1mfSHZUGf8O50Cr0pPZ4zIDjTac).
##' The raw data and the quantification data can also be found in the
##' massIVE repository `MSV000089159`:
##' ftp://massive.ucsd.edu/MSV000089159.
##'
##' @references
##' Andrew Leduc, Gray Huffman, and Nikolai Slavov. 2022. “Droplet
##' Sample Preparation for Single-Cell Proteomics Applied to the Cell
##' Cycle.” bioRxiv. [Link to article](https://doi.org/10.1101/2021.04.24.441211)
##'
##' Andrew Leduc, Gray Huffman, Joshua Cantlon, Saad Khan, and Nikolai
##' Slavov. 2022. “Exploring Functional Protein Covariation across
##' Single Cells Using nPOP.” Genome Biology 23 (1): 261.
##' [Link to article](http://dx.doi.org/10.1186/s13059-022-02817-5)
##'
##' Jason Derks, Andrew Leduc, Georg Wallmann, Gray Huffman, Matthew
##' Willetts, Saad Khan, Harrison Specht, Markus Ralser, Vadim
##' Demichev, and Nikolai Slavov. 2023. “Increasing the Throughput of
##' Sensitive Proteomics by plexDIA.” Nature Biotechnology 41 (1):
##' 50–59. [Link to article](http://dx.doi.org/10.1038/s41587-022-01389-w)
##'
##' @seealso
##' [leduc2022_pSCoPE]
##'
##' @examples
##' \donttest{
##' leduc2022_plexDIA()
##' }
##'
##' @keywords datasets
##'
"leduc2022_plexDIA"

####---- derks2022 ----####

##' Derks et al. 2022 - plexDIA (Nat. Biotechnol.): PDAC vs melanoma
##' cells vs monocytes
##'
##' Single cell proteomics data acquired by the Slavov Lab using the
##' plexDIA protocol. It contains quantitative information from
##' pancreatic ductal acinar cells (PDAC; HPAF-II), melanoma cells
##' (WM989-A6-G3) and monocytes (U-937) at precursor and protein
##' level. The each run acquired 3 samples thanks to mTRAQ
##' multiplexing.
##'
##' @format A [QFeatures] object with 66 assays, each assay being a
##' [SingleCellExperiment] object. The assays either hold the DIA-NN
##' main output report table or the DIA-NN MS1 extracted signal table.
##' The DIA-NN main output report table contains the results of the spectrum
##' identification and quantification. The DIA-NN MS1 extracted
##' signal table contains quantification for all mTRAQ channels if its
##' precursors was identified in at least one of the channels,
##' regardless of whether there is sufficient evidence in those
##' channels at 1% FDR.
##'
##' The data is composed of three datasets
##'
##' 1. **Bulk**: dataset containing bulk (100-cell) data acquired
##'    using a Q-Exactive mass spectrometer. Assays 1-3 contain data
##'    from the DIA-NN main output report; assay 4 is the DIA-NN MS1
##'    extracted signal.
##' 2. **tims**: dataset containing single-cell data acquired using a
##'    timsTOF-SCP mass spectrometer. Assays 5-15  contain data
##'    from the DIA-NN main output report; assay 16 is the DIA-NN MS1
##'    extracted signal.
##' 3. **qe**: dataset containing single-cell data acquired
##'    using a Q-Exactive mass spectrometer. Assays 17-64 contain data
##'    from the DIA-NN main output report; assay 65 is the DIA-NN MS1
##'    extracted signal.
##'
##' The last assay `proteins` contains the processed protein data
##' table generated by the authors.
##'
##' The `colData(derks2022())` contains cell type annotations and
##' batch annotations. The description of the `rowData` fields for the
##' different assays can be found in the
##' [`DIA-NN` documentation](https://github.com/vdemichev/DiaNN#readme).
##'
##' @section Acquisition protocol:
##'
##' The data were acquired using the following setup. More information
##' can be found in the source article (see `References`).
##'
##' - **Cell isolation**: CellenONE cell sorting.
##' - **Sample preparation** performed using the improved SCoPE2
##'   protocol using the CellenONE liquid handling system. nPOP cell
##'   lysis (DMSO) + trypsin digestion + mTRAQ (3plex) labelling and
##'   pooling. A target library was generated as well to
##'   perform prioritized DDA (Huffman et al. 2022) using MaxQuant.Live
##'   (2.0.3).
##' - **Separation**: `bulk` - online nLC (Dionex UltiMate 3000 UHPLC)
##'   with a 25 cm × 75 µm IonOpticks Aurora Series UHPLC column
##'   (AUR2-25075C18A), 200nL/min. `qe` - online nLC (Dionex UltiMate
##'   3000 UHPLC) with a 15 cm × 75 µm IonOpticks Aurora Series UHPLC
##'   column (AUR2-15075C18A), 200nL/min. `tims` - nanoElute liquid
##'   chromatography system (Bruker Daltonics) using a 25 cm × 75 µm,
##'   1.6-µm C18 (AUR2-25075C18A-CSI, IonOpticks).
##' - **Ionization**: ESI.
##' - **Mass spectrometry**: cf article.
##' - **Data analysis**: DIA-NN (1.8.1 beta 16).
##'
##' @section Data collection:
##'
##' The data were collected from a shared Google Drive
##' [folder](https://drive.google.com/drive/folders/1pUC2zgXKtKYn22mlor0lmUDK0frgwL)
##' that is accessible from the SlavovLab website (see `Source` section).
##'
##' For each dataset separately, we combined the sample annotation
##' and the DIANN tables in a [QFeatures] object following the `scp`
##' data structure. We then combined the three datasets in a single
##' `QFeatures` object. We load the proteins table processed by the
##' authors as a [SingleCellExperiment] object and adapted the sample
##' names to match those in the `QFeatures` object. We added the
##' protein data as a new assay and link the precursors to proteins
##' using the `Protein.Group` variable from the `rowData`.
##'
##' @source
##' The data were downloaded from the
##' [Slavov Lab](https://scp.slavovlab.net/Derks_et_al_2022) website.
##' The raw data and the quantification data can also be found in the
##' massIVE repository `MSV000089093`.
##'
##' @references
##' Derks, Jason, Andrew Leduc, Georg Wallmann, R. Gray Huffman,
##' Matthew Willetts, Saad Khan, Harrison Specht, Markus Ralser,
##' Vadim Demichev, and Nikolai Slavov. 2022. "Increasing the
##' Throughput of Sensitive Proteomics by plexDIA." Nature
##' Biotechnology, July.
##' [Link to article](http://dx.doi.org/10.1038/s41587-022-01389-w)
##'
##' @examples
##' \donttest{
##' derks2022()
##' }
##'
##' @keywords datasets
##'
"derks2022"


####---- brunner2022 ----####

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

####---- woo2022_macrophage ----####

##' Woo et al. 2022 (Cell Syst.): LPS-treated macrophages
##'
##' Single-cell data from macrophages subjected to 3 LPS
##' treatments. The data were
##' acquired using the TIFF (transfer identification based on FAIMS
##' filtering) acquisition method. The data contain 155 single cells:
##' 54 control cells (no treatment), 52 cells treated with LPS during
##' 24h and 49 cells treated with LPS during 49h.
##'
##' @format A [QFeatures] object with 5 assays, each assay being a
##' [SingleCellExperiment] object:
##'
##' - `peptides_[intensity or LFQ]`: 2 assays containing peptide
##'   quantities or normalized quantities using the maxLFQ method
##'   as computed by MaxQuant.
##' - `proteins_[intensity or iBAQ or LFQ]`: 3 assays containing
##'   protein quantities or normalized proteins using the iBAQ or
##'   maxLFQ methods as computed by MaxQuant.
##'
##' Sample annotation is stored in `colData(woo_macrophage())`.
##'
##' @section Acquisition protocol:
##'
##' The data were acquired using the following setup. More information
##' can be found in the source article (see `References`).
##'
##' - **Sample isolation**: cultured RAW 264.7 cells treated or not
##'   with 100 ng/ul LPS. The cells were sorted using the Influx II
##'   cell sorter and deposited on a nanoPOTS chip.
##' - **Sample preparation**: cells are lysed using using a DDM+DTT
##'   lysis and reduction buffer. The proteins are alkylated with IAA
##'   and digested with LysC and trypsin. Samples are then acidified
##'   with FA, vacuum dried and stored in freezer until data
##'   acquisition.
##' - **Liquid chromatography**: peptides are loaded using an in-house
##'   autosampler (Williams et al. 2020). The samples are concentrated
##'   through a SPE column (4cm x 100µm i.d. packed with 5µm C18) with
##'   microflow LC pump. The peptides are then eluted from an LC
##'   column (25cm x 50 µm i.d. packed with 1.7µm C18) from a 60 min
##'   gradient (100nL/min).
##' - **Mass spectrometry**: MS/MS was performed on an Orbitrap Fusion
##'   Lumos Tribrid MS with FAIMSpro coupled to a 2.4 kV ESI. FAIMS
##'   setup: 4-CV method (-45, -55, -65, -75 V). MS1 setup: resolution
##'   = 120.000, range = 350-1500 m/z,AGC target of 1E6, accumulation
##'   of 254ms. MS2 setup: 30% HCD, resolution AGC 2E4, accumulation
##'   of 254ms.
##' - **Raw data processing**: preprocessing using Maxquant v1.6.2.10
##'   that use Andromeda search engine (with UniProtKB 2016-21-29).
##'   MBR was enabled.
##'
##' @section Data collection:
##'
##' All data were collected from the MASSIVE repository (accession ID:
##' MSV000085937).
##'
##' The peptide and protein data were extracted from the
##' `peptides_RAW_LPS_scProteomics.txt` or
##' `proteinGroups_RAW_LPS_scProteomics.txt` files, respectively, in
##' the `RAW_LPS_SingleCellProteomics` folders.
##'
##' The tables were split so that intensities, maxLFQ, and iBAQ
##' data are contained in separate tables. Tables are then
##' converted to [SingleCellExperiment] objects. Sample annotations
##' were inferred from the sample names. All data is combined in
##' a [QFeatures] object. [AssayLinks] were stored between peptide
##' assays and their corresponding proteins assays based on the
##' leading razor protein (hence only unique peptides are linked to
##' proteins).
##'
##' The script to reproduce the `QFeatures` object is available at
##' `system.file("scripts", "make-data_woo2022_macrophage.R", package = "scpdata")`
##'
##' @source
##'
##' The peptide and protein data can be downloaded from the MASSIVE
##' repository MSV000085937
##'
##' @references
##'
##' **Source article**: Woo, Jongmin, Geremy C. Clair, Sarah M.
##' Williams, Song Feng, Chia-Feng Tsai, Ronald J. Moore, William B.
##' Chrisler, et al. 2022. “Three-Dimensional Feature Matching
##' Improves Coverage for Single-Cell Proteomics Based on Ion Mobility
##' Filtering.” Cell Systems 13 (5): 426–34.e4.
##' ([link to article](http://dx.doi.org/10.1016/j.cels.2022.02.003)).
##'
##' @examples
##' \donttest{
##' woo2022_macrophage()
##' }
##'
##' @keywords datasets
##'
"woo2022_macrophage"

####---- woo2022_lung ----####

##' Woo et al. 2022 (Cell Syst.): 26 primary human lung cells
##'
##' Single-cell proteomics data from dissociated primary human lung
##' cells. The data were
##' acquired using the TIFF (transfer identification based on FAIMS
##' filtering) acquisition method. The data contain 26 single cells.
##'
##' @format A [QFeatures] object with 5 assays, each assay being a
##' [SingleCellExperiment] object:
##'
##' - `peptides_[intensity or LFQ]`: 2 assays containing peptide
##'   quantities or normalized quantities using the maxLFQ method
##'   as computed by MaxQuant.
##' - `proteins_[intensity or iBAQ or LFQ]`: 3 assays containing
##'   protein quantities or normalized proteins using the iBAQ or
##'   maxLFQ methods as computed by MaxQuant.
##'
##' Sample annotation is stored in `colData(woo_lung())`.
##'
##' @section Acquisition protocol:
##'
##' The data were acquired using the following setup. More information
##' can be found in the source article (see `References`).
##'
##' - **Sample isolation**: primary human lung cells were dissociated
##'   following the protocol in Bandyopadhyay et al., 2018. The cells
##'   were sorted using the Influx II cell sorter and deposited on a
##'   nanoPOTS chip.
##' - **Sample preparation**: cells are lysed using using a DDM+DTT
##'   lysis and reduction buffer. The proteins are alkylated with IAA
##'   and digested with LysC and trypsin. Samples are then acidified
##'   with FA, vacuum dried and stored in freezer until data
##'   acquisition.
##' - **Liquid chromatography**: peptides are loaded using an in-house
##'   autosampler (Williams et al. 2020). The samples are concentrated
##'   through a SPE column (4cm x 100µm i.d. packed with 5µm C18) with
##'   microflow LC pump. The peptides are then eluted from an LC
##'   column (25cm x 50 µm i.d. packed with 1.7µm C18) from a 60 min
##'   gradient (100nL/min).
##' - **Mass spectrometry**: MS/MS was performed on an Orbitrap Fusion
##'   Lumos Tribrid MS with FAIMSpro coupled to a 2.4 kV ESI. FAIMS
##'   setup: 4-CV method (-45, -55, -65, -75 V). MS1 setup: resolution
##'   = 120.000, range = 350-1500 m/z,AGC target of 1E6, accumulation
##'   of 254ms. MS2 setup: 30% HCD, resolution AGC 2E4, accumulation
##'   of 254ms.
##' - **Raw data processing**: preprocessing using Maxquant v1.6.2.10
##'   that use Andromeda search engine (with UniProtKB 2016-21-29).
##'   MBR was enabled.
##'
##' @section Data collection:
##'
##' All data were collected from the MASSIVE repository (accession ID:
##' MSV000085937).
##'
##' The peptide and protein data were extracted from the
##' `peptides_nondepleted_Lung_scProteomics.txt` or
##' `proteinGroups_nondepleted_Lung_scProteomics.txt` files,
##' respectively, in the `NonDepleted_Lung_SingleCellProteomics`
##' folders.
##'
##' The tables were split so that intensities, maxLFQ, and iBAQ
##' data are contained in separate tables. Tables are then
##' converted to [SingleCellExperiment] objects. Sample annotations
##' were inferred from the sample names. All data is combined in
##' a [QFeatures] object. [AssayLinks] were stored between peptide
##' assays and their corresponding proteins assays based on the
##' leading razor protein (hence only unique peptides are linked to
##' proteins).
##'
##' The script to reproduce the `QFeatures` object is available at
##' `system.file("scripts", "make-data_woo2022_lung.R", package = "scpdata")`
##'
##' @source
##'
##' The peptide and protein data can be downloaded from the MASSIVE
##' repository MSV000085937
##'
##' @references
##'
##' **Source article**: Woo, Jongmin, Geremy C. Clair, Sarah M.
##' Williams, Song Feng, Chia-Feng Tsai, Ronald J. Moore, William B.
##' Chrisler, et al. 2022. “Three-Dimensional Feature Matching
##' Improves Coverage for Single-Cell Proteomics Based on Ion Mobility
##' Filtering.” Cell Systems 13 (5): 426–34.e4.
##' ([link to article](http://dx.doi.org/10.1016/j.cels.2022.02.003)).
##'
##' @examples
##' \donttest{
##' woo2022_lung()
##' }
##'
##' @keywords datasets
##'
"woo2022_lung"

####---- gregoire2023_mixCTRL ----####

##' Grégoire et al. 2023 - mixCTRL (arXiv): benchmark using
##' monocytes/macrophages
##'
##' Single cell proteomics data acquired using the SCoPE2 protocol.
##' The dataset contains two monocytes cell lines (THP1 and U937) as
##' well as controled mixtures of both and macrophage-like cells
##' produced upon PMA treatment. It contains quantitative information
##' at PSM, peptide and protein levels. Data was acquired using Lumos
##' Orbitrap (mainly) and timsTOF SCP mass spectrometers.
##'
##' @format A [QFeatures] object with 119 assays, each assay being a
##' [SingleCellExperiment] object:
##'
##' - Assays 1-42: PSM data acquired with a TMT-16plex protocol, hence
##'   those assays contain 16 columns. Columns hold quantitative
##'   information from single-cell channels, carrier channels,
##'   blank (negative control) channels and unused channels.
##' - Assays 43-84: peptide data resulting from the PSM to peptide
##'   aggregation of the 42 PSM assays.
##' - Assays 85-91: peptide data for each of the 7 acquisition
##'   batches. Peptide data were joined based on their respective
##'   acquisition batches.
##' - Assays 92-98: normalised peptide data.
##' - Assays 99-105: normalised and log-transformed peptide data.
##' - Assays 106-112: protein data for each of the 7 acquisition
##'   batches. Normalised and log-transformed peptide data were
##'   agreggated to protein.
##' - Assays 113-119: Batch corrected protein data. Normalised and
##'   log-transformed protein data were batch corrected to remove
##'   technical variability induced by runs and channels.
##'
##' All the data has been filtered to keep high quality features and
##' samples.
##'
##' The `colData(gregoire2023_mixCTRL())` contains cell type annotation and
##' batch annotation that are common to all assays. The description of
##' the `rowData` fields for the PSM data can be found in the
##' [`sage` documentation](https://sage-docs.vercel.app/docs/results/search).
##'
##' @section Acquisition protocol:
##'
##' The data were acquired using the following setup. More information
##' can be found in the source article (see *References*).
##'
##' - **Cell isolation**: BD FACSAria III cell sorting.
##' - **Sample preparation** performed using the SCoPE2 protocol: mPOP
##'   cell lysis + trypsin digestion + TMT-16plex labeling and
##'   pooling.
##' - **Separation**: online nLC (Ultimate 3000 LC System or Vanquish
##'   Neo UHPLC System) with a BioZen Peptide Polar C18 250 x 0.0075mm
##'   column.
##' - **Mass spectrometry**:  Orbitrap Fusion Lumos Tribrid (MS1
##'   resolution = 70,000; MS2 accumulation time = 120ms; MS2
##'   resolution = 70,000) and timsTOF SCP.
##' - **Data preprocessing**: Sage.
##'
##' @section Data collection:
##'
##' The PSM data were collected from a Zenodo archive (see `Source`
##' section). The folder contains the following files of interest:
##'
##' - `results.sage.cbio.tsv`: the sage identification output file for
##'   batches acquired on the Lumos MS.
##' - `results.sage.giga.tsv`: the sage identification output file for
##'   batches acquired on the timsTOF SCP MS.
##' - `quant.cbio.tsv`: the sage quantification output file for
##'   batches acquired on the Lumos MS.
##' - `quant.giga.tsv`: the sage quantification output file for
##'   batches acquired on the timsTOF SCP MS.
##' - `sampleAnnotation_batch.csv`: sample annotation for each
##'   acquisition batch. There are in total 8 different annotation
##'   files.
##'
##' We combined the sample annotations in a single table. We also
##' combined `cbio` and `giga` tables together and merged resulting
##' identification and quantification tables. Both annotation and
##' features tables are then combined in a single [QFeatures] object
##' using the [scp::readSCP()] function.
##'
##' The [QFeatures] object was processed as described in the author's
##' manuscript (see `source`). Note that the imputed assays were used
##' in the paper for illustrative purposes only and have not been
##' reproduced here.
##'
##' @source
##' The data were downloaded from the [Zenodo
##' repository](https://zenodo.org/records/8417228).  The raw data and
##' the quantification data can also be found in the ProteomeXchange
##' Consortium via the [PRIDE partner
##' repository](https://www.ebi.ac.uk/pride/archive/projects/PXD046211),
##' project `PXD046211`.
##'
##' @references
##' Samuel Grégoire, Christophe Vanderaa, Sébastien Pyr dit Ruys,
##' Gabriel Mazzucchelli, Christopher Kune, Didier Vertommen and
##' Laurent Gatto. 2023. *Standardised workflow for mass spectrometry-
##' based single-cell proteomics data processing and analysis using
##' the scp package.*
##' arXiv. DOI:[10.48550/arXiv.2310.13598](https://doi.org/10.48550/arXiv.2310.13598)
##'
##' @examples
##' \donttest{
##' gregoire2023_mixCTRL()
##' }
##'
##' @keywords datasets
##'
"gregoire2023_mixCTRL"

####---- khan2023 ----####


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
##' - `proteins_imputed`: protein data containing quantitative data for 4096
##'   proteins and 421 single-cells with k-nearest neighbors (KNN) imputation.
##' - `proteins_unimputed`: protein data containing quantitative data for 4096
##'   proteins and 421 single-cells without imputation.
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
##' (`EpiToMesen.TGFB.nPoP_trial1_ProtByCellMatrix_NSThreshDART_medIntCrNorm_imputedNotBC.csv`).
##' The data were formatted to a [SingleCellExperiment] object and the sample
##' metadata were matched to the column names (mapping is retrieved
##' after running the SCoPE2 R script, `EMTTGFB_singleCellProcessing.R`) and
##' stored in the `colData`. The object is then added to the [QFeatures] object
##' and the rows of the peptide data are linked to the rows of the protein data
##' based on the protein sequence information through an `AssayLink` object.
##'
##' The unimputed protein data were taken from the same google drive folder
##' (`EpiToMesen.TGFB.nPoP_trial1_ProtByCellMatrix_NSThreshDART_medIntCrNorm_unimputed.csv`).
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
##' transition." bioRxiv.
##' ([link to article](https://doi.org/10.1101/2023.12.21.572913)).
##'
##' @examples
##' \donttest{
##' khan2023()
##' }
##'
##' @keywords datasets
##'
"khan2023"


####---- guise2024 ----####


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

####---- petrosius2023_mES ----####

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
##'   using a 130 μm sorting chip. Cells were sorted at single-cell resolution, 
##'   into a 384-well Eppendorf LoBind PCR plate (Eppendorf AG) containing 1 μL 
##'   of lysis buffer.
##' - **Sample preparation**: Single-cell protein lysates were digested with 
##'   2 ng of Trypsin (Sigma cat. Nr. T6567) supplied in 1 μL of digestion 
##'   buffer (100mM TEAB pH 8.5, 1:5000 (v/v) benzonase (Sigma cat. Nr. E1014)).
##'   The digestion was carried out overnight at 37 °C, and subsequently 
##'   acidified by the addition of 1 μL 1% (v/v) trifluoroacetic acid (TFA). 
##'   All liquid dispensing was done using an I-DOT One instrument (Dispendix).
##' - **Liquid chromatography**: The Evosep one liquid chromatography system was 
##'   used for DIA isolation window survey and HRMS1-DIA experiments.The standard
##'   31 min or 58min pre-defined Whisper gradients were used, where peptide 
##'   elution is carried out with 100 nl/min flow rate. A 15 cm × 75 μm 
##'   ID column (PepSep) with 1.9 μm C18 beads (Dr. Maisch, Germany) and a 10 
##'   μm ID silica electrospray emitter (PepSep) was used. Both LC systems were 
##'   coupled online to an orbitrap Eclipse TribridMass Spectrometer 
##'   (ThermoFisher Scientific) via an EasySpray ion source connected to a 
##'   FAIMSPro device.
##' - **Mass spectrometry**: The mass spectrometer was operated in positive 
##'   mode with the FAIMSPro interface compensation voltage set to −45 V.
##'   MS1 scans were carried out at 120,000 resolution with an automatic gain
##'   control (AGC) of 300% and maximum injection time set to auto. For the DIA 
##'   isolation window survey a scan range of 500–900 was used and 400–1000 
##'   rest of the experiments. Higher energy collisional dissociation (HCD) was 
##'   used for precursor fragmentation with a normalized collision energy (NCE) 
##'   of 33% and MS2 scan AGC target was set to 1000%. 
##' - **Raw data processing**: The mESC raw data files were processed with 
##'   Spectronaut 17 and protein abundance tables exported and analyzed further 
##'   with python. 
##'
##' @section Data collection:
##'
##' The data were provided by the Author and is accessible at the [Dataverse]
##' (https://dataverse.uclouvain.be/dataset.xhtml?persistentId=doi:10.14428/DVN/EMAVLT)
##' The folder ('20240205_111248_mESC_SNEcombine_m15-m2i/') contains the 
##' following files of interest:
##'
##' - `20240205_111251_PEPQuant (Normal).tsv`: the PSM level data
##' - `20240205_111251_Peptide Quant (Normal).tsv`: the peptide level data
##' - `20240205_111251_PGQuant (Normal).tsv`: the protein level data
##'
##' The metadata were downloaded from the [Zenodo
##' repository] (https://zenodo.org/records/8146605).
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
##' @source
##' The peptide and protein data can be downloaded from the [Dataverse]
##' (https://dataverse.uclouvain.be/dataset.xhtml?persistentId=doi:10.14428/DVN/EMAVLT)
##' The raw data and the quantification data can also be found in the
##' MassIVE repository `MSV000092429`:
##' ftp://MSV000092429@massive.ucsd.edu/.
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

####---- petrosius2023_AML ----####

##' Petrosius et al. 2023 (bioRxiv): AML hierarchy on Astral.
##'
##' Single cell proteomics data from FACS sorted cells from the
##' OCI-AML8227 model. The dataset contains leukemic stem cells (LSC;
##' CD34+, CD38-), progenitor cells (CD34+, CD38+), CD38+ blasts
##' (CD34-, CD38+) and CD38- blasts (CD34-, CD38-). It contains
##' quantitative information at PSM, peptide and protein levels. Data
##' was acquired using an Orbitrap Astral mass spectrometer. Direct DIA
##' analysis was performed with Spectronaut version 17.
##'
##' @format A [QFeatures] object with 217 assays, each assay being a
##' [SingleCellExperiment] object:
##'
##' - Assays 1-215: PSM data from the Spectronaut PEPQuant file with
##'   LFQ quantities from the FG.MS1Quantity column.
##' - `peptides`: Peptide data resulting from the PSM to peptide
##'   aggregation the 215 PSM assays. Resulting peptide assays were
##'   joined into a single assay.
##' - `proteins`: Protein data from the Spectronaut PGQuant file with
##'   LFQ quantities from the PG.Quantity column.
##'
##' The `colData(petrosius2023_AML())` contains cell type annotation, batch
##' annotation and FACS data. The description of the `rowData` fields
##' can be found in the [`Spectronaut` user manual](https://biognosys.com/content/uploads/2023/03/Spectronaut-17_UserManual.pdf).
##'
##' @section Acquisition protocol:
##'
##' The data were acquired using the following setup. More information
##' can be found in the source article (see *References*).
##'
##' - **Cell isolation**: Cell sorting was done on a FACS Aria III or
##'   Aria II instrument, controlled by the DIVA software package and
##'   operated with a 100 μm nozzle. Cells were sorted at single-cell
##'   resolution, into a 384-well Eppendorf LoBind PCR plate containing
##'   1 μL of lysis buffer.
##' - **Sample preparation** Single-cell protein lysates were digested
##'   overnight at 37°C with 2 ng of Trypsin supplied in 1 μL of
##'   digestion buffer. Digestion was stopped by the addition of 1 μL
##'   1% (v/v) trifluoroacetic acid (TFA). All liquid dispensing was
##'   done using an I-DOT One instrument.
##' - **Liquid chromatography**: Chromatographic separation of peptides
##'   was conducted on a vanquish Neo UHPLC system connected to a 50 cm
##'   uPAC Neo Low-load and an EASY-spray. Autosampler and injection
##'   valves were configured to perform direct injections from a 384
##'   well plate using a 25 uL injection loop on 11.8 min gradients.
##' - **Mass spectrometry**: Acquisition was conducted with an Orbitrap
##'   Astral mass spectrometer operated in positive mode with the
##'   FAIMSPro interface compensation voltage set to −45 V.
##'   MS1 scans were acquired with the Orbitrap at a resolution of
##'   120,000 and a scan range of 400 to 900 m/z with normalized
##'   automatic gain control (AGC) target of 300 % and maximum
##'   injection time of 246 ms. Data independent acquisition of MS2
##'   spectra was performed in the Astral using loop control set to 0.7
##'   seconds per cycle with varying isolation window widths and
##'   injection times. Fragmentation of precursor ions was performed
##'   using higher energy collisional dissociation (HCD) using a
##'   normalized collision energy (NCE) of 25 %. AGC target was set to
##'   800 %.
##' - **Raw data processing**: Raw files were processed using
##'   Spectronaut version 17. Direct DIA analysis was performed in
##'   pipeline mode. Pulsar searches were performed without fixed
##'   modifications. N-terminal acetylation and methionine oxidation
##'   were set as variable modifications. Quantification level was set
##'   to MS1 and the quantity type set to area under the curve.
##'
##' @section Data collection:
##'
##' The data were provided by the authors and is accessible at the [Dataverse]
##' (https://dataverse.uclouvain.be/dataset.xhtml?persistentId=doi:10.14428/DVN/EMAVLT)
##' The dataset ('Astral AML single-cell data from Petrosius et al. 2023 preprint')
##' contains the following files of interest:
##'
##' - `20240201_130747_PEPQuant (Normal).tsv`: the PSM level data
##' - `20240201_130747_PGQuant (Normal).tsv`: the protein level data
##' - `index_map.csv`: FACS data.
##' - `msRuns_overview.csv`: Sample annotations.
##'
##' We added the FACS data to the sample annotations in a single table.
##' Both annotations and PSM features tables are then combined in a
##' single [QFeatures] object using the [scp::readSCP()] function.
##'
##' The peptide data were obtained by aggregation of the PSM data to
##' the peptide level. All of the resulting peptides assays were joined
##' into a single assays. Individual peptides assays were discarded.
##'
##' The protein data were formatted from the `20240201_130747_PGQuant (Normal).tsv`
##' to a [SingleCellExperiment] object and the sample metadata were
##' matched to the column names and stored in the  `colData`. The
##' object is then added to the [QFeatures] object and the rows of the
##' peptide data are linked to the rows of the protein data based on
##' the protein sequence information through an `AssayLink` object.
##'
##' Note that the [QFeatures] object has not been further processed and
##' has therefore not been normalized, log-transformed or
##' batch-corrected.
##'
##' @source
##' The PSM data, protein data and sample annotations can be downloaded
##' from the dataset 'Astral AML single-cell data from Petrosius et al. 2023 preprint'
##' in the [Dataverse](https://dataverse.uclouvain.be/dataset.xhtml?persistentId=doi:10.14428/DVN/EMAVLT).
##'
##' @references
##' Valdemaras Petrosius, Pedro Aragon-Fernandez, Tabiwang N. Arrey,
##' Nil Üresin, Benjamin Furtwängler, Hamish Stewart, Eduard Denisov,
##' Johannes Petzoldt, Amelia C. Peterson, Christian Hock, Eugen Damoc,
##' Alexander Makarov, Vlad Zabrouskov, Bo T. Porse and Erwin M. Schoof.
##' 2023. "Evaluating the capabilities of the Astral mass analyzer for
##' single-cell proteomics." biorxiv.
##' https://doi.org/10.1101/2023.06.06.543943
##' DOI:[10.1101/2023.06.06.543943](https://doi.org/10.1101/2023.06.06.543943)
##'
##' @examples
##' \donttest{
##' petrosius2023_AML()
##' }
##'
##' @keywords datasets
##'
"petrosius2023_AML"
