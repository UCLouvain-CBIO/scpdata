

####---- SPECHT ET AL. 2019 ----####


##' Specht et al. 2019: macrophages vs monocytes (version 2)
##'
##' Single cell proteomics data acquired by the Slavov Lab. This is the version 
##' 2 of the data released in December 2019. It contains quantitative 
##' information of macrophages and monocytes at PSM, peptide and protein level. 
##' 
##' @format A `Features` object with 179 assays, each assay being a 
##' `SingleCellExperiment` object: 
##' 
##' - Assay 1-63: PSM data for SCoPE2 sets acquired with a TMT-11 
##'   multiplexing protocole, hence those assays contain 11 columns. Columns 
##'   hold quantitative information from single-cell channels, carrier channels, 
##'   reference channels, empty (blank) channels and unused channels.
##' - Assay 64-177: PSM data for SCoPE2 sets acquired with a TMT-16
##'   multiplexing protocole, hence those assays contain 16 columns. Columns 
##'   hold quantitative information from single-cell channels, carrier channels, 
##'   reference channels, empty (blank) channels and unused channels.
##' - `peptides`: peptide data containing quantitative data for 9208 
##'   peptides and 1018 single-cells. Cell type annotation and batch annotation
##'   are stored in `colData(specht2019v2[[178]]`.
##' - `proteins`: protein data containing quantitative data for 2772 
##'   proteins and 1018 single-cells. Cell type annotation and batch annotation
##'   are stored in `colData(specht2019v2[[179]]`.
##' }
##' The `colData(specht2019v2)` contains cell type annotation and batch 
##' annotation that are common to all assays.
##' 
##' See `Details`` for information about data collection.
##'
##' @details 
##' 
##' **Acquisition protocole**
##' 
##' The data was acquired using the following setup. More information can be 
##' found in the source article (see `References`).
##' 
##' - **Cell isolation**: flow cytometry (BD FACSAria I).
##' - **Sample preparation** performed using the SCoPE2 protcole. mPOP cell 
##'   lysis + trypsin digestion + TMT 11plex or 16plex labeling and pooling.
##' - **Separation**: online nLC (DionexUltiMate 3000 UHPLC with a 25cm x 75μm 
##'   IonOpticksAurora Series UHPLC column; 200nL/min).
##' - **Ionization**: ESI (2,200V).
##' - **Mass spectrometry**: Thermo Scientific Q-Exactive (MS1 resolution = 
##'   70,000; MS2 accumulation time = 300ms; MS2 resolution = 70,000).
##' - **Data analysis**: DART-ID + MaxQuant (1.6.2.3).
##' 
##' **Data collection**
##' 
##' The PSM data were collected from a shared Goolge Drive folder that 
##' is accessible from the SlavovLab website (see `Source` section). The folder 
##' contains the following 
##' files of interest: 
##' 
##' - `ev_updated.txt`: the MaxQuant/DART-ID output file
##' - `annotation_fp60-97.csv`: sample annotation
##' - `batch_fp60-97.csv`: batch annotation
##' 
##' We combined the the sample annotation and the batch annotation in a single 
##' table. We also formated the quantification table so that columns match with 
##' those of the annotation and filter only for single-cell runs. Both table 
##' are then combined in a single `Features` object, where the quantitative data
##' are split with respect to batch. 
##'  
##' The peptide data were taken from the Slavov lab directly (`Peptides-raw.csv`). 
##' It is provided as a spreadsheet. The data were formated to a 
##' `SingleCellExperiment` object and the sample metadata were matched to the 
##' column names (mapping is retrieved after running the SCoPE2 R script) and 
##' stored in the `colData`. The object is then added to the `Features` object 
##' (containing the PSM assays) and the rows of the peptide data are linked to 
##' the rows of the PSM data based on the peptide sequence information through 
##' an `AssayLink` object. 
##' 
##' The protein data (`Proteins-processed.csv`) is formated similarly to the 
##' peptide data, and the rows of the proteins were mapped onto the rows of the 
##' peptide data based on the protein sequence information.
##'  
##' @source 
##' The data were downloaded from the 
##' [Slavov Lab](https://scope2.slavovlab.net/docs/data) website via a
##' shared Google Drive 
##' [folder](https://drive.google.com/drive/folders/1VzBfmNxziRYqayx3SP-cOe2gu129Obgx). 
##' The raw data and the quantification data can also be found in the massIVE 
##' repository 
##' [MSV000083945](ftp://massive.ucsd.edu/MSV000083945).
##' 
##' @references Specht, Harrison, Edward Emmott, Aleksandra A. Petelski, R. Gray 
##' Huffman, David H. Perlman, Marco Serra, Peter Kharchenko, Antonius Koller, 
##' and Nikolai Slavov. 2019. “Single-Cell Mass-Spectrometry Quantifies the 
##' Emergence of Macrophage Heterogeneity.” bioRxiv. 
##' ([link to article](https://doi.org/10.1101/665307)).
##' 
##' @docType data
##'
##' @keywords datasets
##' 
"specht2019v2"


####---- DOU ET AL. 2019 ----####


##' Dou et al. 2019: HeLa digests 
##'  
##' Single-cell proteomics using nanoPOTS combined with TMT isobaric labeling.
##' It contains quantitative information at PSM and protein level.
##' The samples are commercial Hela digest diluted to single cell amounts 
##' (0.2ng). The boosting wells contain the same digest but at higher dose (10 
##' ng).
##' 
##' @format A `Features` object with 3 assays, each assay being a 
##' `SingleCellExperiment` object: 
##' 
##' - `Hela_run_1`: PSM data with 10 columns corresponding to the TMT 10plex 
##'   channels. Columns hold quantitative information for HeLa digestion samples
##'   (either 0, 0.2 or 10ng). This is the data for run 1. 
##' - `Hela_run_1`: PSM data with 10 columns corresponding to the TMT 10plex 
##'   channels. Columns hold quantitative information for HeLa digestion samples
##'   (either 0, 0.2 or 10ng). This is the data for run 2.
##' - `proteins`: protein data containing quantitative data for 1641 proteins 
##'   and 20 samples (run 1 and run 2 combined). 
##' }
##' 
##' Sample annotation is stored in `colData(dou2019_hela)`.
##' 
##' See `Details` for information about data collection.
##' 
##' @details 
##' 
##' **Acquisition protocole**
##' 
##' The data was acquired using the following setup. More information can be 
##' found in the source article (see `References`).
##' 
##' - **Cell isolation**: commercially available HeLa protein digest (Thermo 
##'   Scientific).
##' - **Sample preparation** performed using the nanoPOTs device. Protein 
##'   extraction (DMM + TCEAP) + alkylation (IAA) + Lys-C digestion + trypsin
##'   digestion + TMT 10plex labeling and pooling.
##' - **Separation**: nanoLC (Dionex UltiMate with an in-house packed 50cm x 
##'   30μm LC columns; 50nL/min)
##' - **Ionization**: ESI (2,000V)
##' - **Mass spectrometry**: Thermo Fisher Orbitrap Fusion Lumos Tribrid (MS1
##'   accumulation time = 50ms; MS1 resolution = 120,000; MS1 AGC = 1E6; MS2 
##'   accumulation time = 246ms; MS2 resolution = 60,000; MS2 AGC = 1E5)
##' - **Data analysis**: MS-GF+ + MASIC (v3.0.7111) + RomicsProcessor (custom R
##'   package)
##' 
##' **Data collection**
##' 
##' The PSM data were collected from the MassIVE repository MSV000084110 (see 
##' `Source` section). The files downloaded files are:
##' 
##' - `Hela_run_*_msgfplus.mzid`: the MS-GF+ identification result files
##' - `Hela_run_*_ReporterIons.txt`: the MASIC quantification result files
##' 
##' For each batch, the quantification and identification data were combined 
##' based on the scan number (common to both data sets). The combined datasets 
##' for the different runs were then concatenated feature-wise. The sample 
##' annotation table was manually created based on the available information 
##' provided in the article. The data were then converted to a `QFeatures` 
##' object using the `scp::readSCP` function. 
##' 
##' The protein data was download from `Supporting information` section from the 
##' publisher's website (see `Sources`). The data is supplied as an Excel file 
##' `ac9b03349_si_003.xlsx`. The file contains 7 sheets from which we only took 
##' the sheet 6 (named `5 - Run 1 and 2 raw data``) with the combined protein 
##' data for the two runs. We converted the data to a `SingleCellExperiment` 
##' object and added the object as a new assay in the `QFeatures` dataset 
##' (containing the PSM data). Links between the proteins and the corresponding
##' PSM were created. 
##' 
##' @source 
##' The PSM data can be downloaded from the massIVE repository MSV000084110.
##' 
##' The protein data can be downloaded from the 
##' [ACS Publications](https://pubs.acs.org/doi/10.1021/acs.analchem.9b03349)
##' website (Supporting information section).
##' 
##' @references 
##' Dou, Maowei, Geremy Clair, Chia-Feng Tsai, Kerui Xu, William B. 
##' Chrisler, Ryan L. Sontag, Rui Zhao, et al. 2019. “High-Throughput Single Cell 
##' Proteomics Enabled by Multiplex Isobaric Labeling in a Nanodroplet Sample 
##' Preparation Platform.” Analytical Chemistry, September 
##' ([link to article](https://doi.org/10.1021/acs.analchem.9b03349)).
##' 
##' @seealso 
##' [dou2019_mouse], [dou2019_boosting]
##' 
##' @docType data
##' 
##' @keywords datasets
##' 
##' 
"dou2019_hela"


##' Dou et al. 2019: single cells from cultured murine cell lines
##'  
##' Single-cell proteomics using nanoPOTS combined with TMT isobaric labeling. 
##' It contains quantitative information at PSM and protein level.
##' The cell types are either "Raw" (macrophage cells), "C10" (epihelial cells), 
##' or "SVEC" (endothelial cells). Out of the 132 wells, 72 contain single cells, 
##' corresponding to 24 C10 cells, 24 RAW cells, and 24 SVEC. The other wells are 
##' either boosting channels (12), empty channels (36) or reference channels 
##' (12). Boosting and reference channels are balanced (1:1:1) mixes of C10, 
##' SVEC, and RAW samples at 5 ng and 0.2 ng, respectively. The different cell 
##' types where evenly distributed across 4 nanoPOTS chips. Samples were 
##' 11-plexed with TMT labeling.
##' 
##' @format A `QFeatures` object with 13 assays, each assay being a 
##' `SingleCellExperiment` object: 
##' 
##' - `Single_Cell_Chip_X_Y`: PSM data with 11 columns corresponding to the TMT 
##'   channels (see `Notes`). The `X` indicates the chip number (from 1 to 4) 
##'   and `Y` indicates the row name on the chip (from A to C). 
##' - `proteins`: protein data containing quantitative data for 2331 proteins 
##'   and 132 samples (all runs combined). 
##' }
##' 
##' Sample annotation is stored in `colData(dou2019_mouse)`.
##' 
##' See `Details` for information about data collection.
##' 
##' @details 
##' 
##' **Acquisition protocole**
##' 
##' The data was acquired using the following setup. More information can be 
##' found in the source article (see `References`).
##' 
##' - **Cell isolation**: single-cells from the three murine cell lines were
##'   isolated using FACS (BD Influx II cell sorter ). 
##' - **Sample preparation** performed using the nanoPOTs device. Protein 
##'   extraction (DMM + TCEAP) + alkylation (IAA) + Lys-C digestion + trypsin
##'   digestion + TMT 10plex labeling and pooling.
##' - **Separation**: nanoLC (Dionex UltiMate with an in-house packed 50cm x 
##'   30μm LC columns; 50nL/min)
##' - **Ionization**: ESI (2,000V)
##' - **Mass spectrometry**: Thermo Fisher Orbitrap Fusion Lumos Tribrid (MS1
##'   accumulation time = 50ms; MS1 resolution = 120,000; MS1 AGC = 1E6; MS2 
##'   accumulation time = 246ms; MS2 resolution = 60,000; MS2 AGC = 1E5)
##' - **Data analysis**: MS-GF+ + MASIC (v3.0.7111) + RomicsProcessor (custom R
##'   package)
##' 
##' **Data collection**
##' 
##' The PSM data were collected from the MassIVE repository MSV000084110 (see 
##' `Source` section). The files downloaded files are:
##' 
##' 
##' - `Single_Cell_Chip_*_*_msgfplus.mzid`: the MS-GF+ identification result 
##'   files.
##' - `Single_Cell_Chip_*_*_ReporterIons.txt`: the MASIC quantification result 
##'   files.
##' 
##' For each batch, the quantification and identification data were combined 
##' based on the scan number (common to both data sets). The combined datasets 
##' for the different runs were then concatenated feature-wise. The sample 
##' annotation table was manually created based on the available information 
##' provided in the article. The data were then converted to a `QFeatures` 
##' object using the `scp::readSCP` function. 
##' 
##' The protein data was download from `Supporting information` section from the 
##' publisher's website (see `Sources`). The data is supplied as an Excel file 
##' `ac9b03349_si_005.xlsx`. The file contains 7 sheets from which we only took 
##' the 2nd (named `01 - Raw sc protein data`) with the combined protein data 
##' for the 12 runs. We converted the data to a `SingleCellExperiment` object 
##' and added the object as a new assay in the `QFeatures` dataset (containing 
##' the PSM data). Links between the proteins and the corresponding PSM were 
##' created. 
##' 
##' @note Although a TMT 10plex labeling is reported in the article, the PSM 
##' data contained 11 channels for each run. Those 11th channel contain mostly 
##' missing data and are hence assumed to be empty channels. 
##' 
##' @source 
##' The PSM data can be downloaded from the massIVE repository MSV000084110.
##' 
##' The protein data can be downloaded from the 
##' [ACS Publications](https://pubs.acs.org/doi/10.1021/acs.analchem.9b03349)
##' website (Supporting information section).
##' 
##' @references 
##' Dou, Maowei, Geremy Clair, Chia-Feng Tsai, Kerui Xu, William B. 
##' Chrisler, Ryan L. Sontag, Rui Zhao, et al. 2019. “High-Throughput Single Cell 
##' Proteomics Enabled by Multiplex Isobaric Labeling in a Nanodroplet Sample 
##' Preparation Platform.” Analytical Chemistry, September 
##' ([link to article](https://doi.org/10.1021/acs.analchem.9b03349)).
##' 
##' @seealso 
##' [dou2019_hela], [dou2019_boosting]
##' 
##' @docType data
##' 
##' @keywords datasets
##' 
##' 
"dou2019_mouse"


##' Dou et al. 2019: testing boosting ratios
##'
##' Single-cell proteomics using nanoPOTS combined with TMT isobaric labeling. 
##' It contains quantitative information at PSM and protein level.
##' The cell types are either "Raw" (macrophage cells), "C10" 
##' (epihelial cells), or "SVEC" (endothelial cells). Each cell is replicated 2 
##' or 3x. Each cell type was run using 3 levels of boosting: 0 ng (no boosting), 
##' 5 ng or 50 ng. When boosting was applied, 1 reference well and 1 boosting 
##' well were added, otherwise 1 empty well was added. Each boosting setting 
##' (0ng, 5ng, 50ng) was run in duplicate.
##' 
##' @format A `QFeatures` object with 7 assays, each assay being a 
##' `SingleCellExperiment` object: 
##' 
##' - `Boosting_X_run_Y`: PSM data with 10 columns corresponding to the TMT 
##'   10plex channels. The `X` indicates the boosting amount (0ng, 5ng or 50ng) 
##'   and `Y` indicates the run number (1 or 2). 
##' - `proteins`: protein data containing quantitative data for 1436 proteins 
##'   and 60 samples (all runs combined). 
##' }
##' 
##' Sample annotation is stored in `colData(dou2019_boosting)`.
##' 
##' See `Details` for information about data collection.
##' 
##' @details 
##' 
##' **Acquisition protocole**
##' 
##' The data was acquired using the following setup. More information can be 
##' found in the source article (see `References`).
##' 
##' - **Cell isolation**: single-cells from the three murine cell lines were
##'   isolated using FACS (BD Influx II cell sorter ). Boosting sample were 
##'   prepared (presumably in bulk) from 1:1:1 mix of the three cell lines.
##' - **Sample preparation** performed using the nanoPOTs device. Protein 
##'   extraction (DMM + TCEAP) + alkylation (IAA) + Lys-C digestion + trypsin
##'   digestion + TMT 10plex labeling and pooling.
##' - **Separation**: nanoLC (Dionex UltiMate with an in-house packed 50cm x 
##'   30μm LC columns; 50nL/min)
##' - **Ionization**: ESI (2,000V)
##' - **Mass spectrometry**: Thermo Fisher Orbitrap Fusion Lumos Tribrid (MS1
##'   accumulation time = 50ms; MS1 resolution = 120,000; MS1 AGC = 1E6; MS2 
##'   accumulation time = 246ms; MS2 resolution = 60,000; MS2 AGC = 1E5)
##' - **Data analysis**: MS-GF+ + MASIC (v3.0.7111) + RomicsProcessor (custom R
##'   package)
##' 
##' **Data collection**
##' 
##' The PSM data were collected from the MassIVE repository MSV000084110 (see 
##' `Source` section). The files downloaded files are:
##' 
##' - `Boosting_*ng_run_*_msgfplus.mzid`: the MS-GF+ identification result 
##'   files.
##' - `Boosting_*ng_run_*_ReporterIons.txt`: the MASIC quantification result 
##'   files.
##' 
##' For each batch, the quantification and identification data were combined 
##' based on the scan number (common to both data sets). The combined datasets 
##' for the different runs were then concatenated feature-wise. The sample 
##' annotation table was manually created based on the available information 
##' provided in the article. The data were then converted to a `QFeatures` 
##' object using the `scp::readSCP` function. 
##' 
##' The protein data was download from `Supporting information` section from the 
##' publisher's website (see `Sources`). The data is supplied as an Excel file 
##' `ac9b03349_si_004.xlsx`. The file contains 7 sheets from which we took the 
##' 2nd, 4th and 6th sheets (named `01 - No Boost raw data`, 
##' `03 - 5ng boost raw data`, `05 - 50ng boost raw data`, respectively). The 
##' sheets contain the combined protein data for the duplicate runs given the 
##' boosting amount. We joined the data for all boosting ration based on the 
##' protein name and converted the data to a `SingleCellExperiment` object. We
##' then added the object as a new assay in the `QFeatures` dataset (containing 
##' the PSM data). Links between the proteins and the corresponding PSM were 
##' created. 
##' 
##' @source 
##' The PSM data can be downloaded from the massIVE repository MSV000084110.
##' 
##' The protein data can be downloaded from the 
##' [ACS Publications](https://pubs.acs.org/doi/10.1021/acs.analchem.9b03349)
##' website (Supporting information section).
##' 
##' @references 
##' Dou, Maowei, Geremy Clair, Chia-Feng Tsai, Kerui Xu, William B. 
##' Chrisler, Ryan L. Sontag, Rui Zhao, et al. 2019. “High-Throughput Single Cell 
##' Proteomics Enabled by Multiplex Isobaric Labeling in a Nanodroplet Sample 
##' Preparation Platform.” Analytical Chemistry, September 
##' ([link to article](https://doi.org/10.1021/acs.analchem.9b03349)).
##' 
##' @seealso 
##' [dou2019_hela], [dou2019_mouse]
##' 
##' @docType data
##' 
##' @keywords datasets
##' 
##' 
"dou2019_boosting"


####---- ZHU ET AL. 2018 (MCP) ----####


##' Zhu et al. 2018 (Mol. Cel. Prot.): rat brain laser dissections
##'
##' Near single-cell proteomics data of laser cuptured micro-dissection samples.
##' The samples are 24 brain section from rat pups (day 17). The slices are 12
##' µm thick squares of either 50, 100, or 200 µm width. 5 samples were dissected 
##' from the corpus callum (`CC`), 4 samples were dissected from the 
##' corpus collosum (`CP`), 13 samples were extracted from the cerebral 
##' cortex (`CTX`), and 2 samples are labeled as (`Mix``).
##' 
##' @format A `QFeatures` object with 1 assay, `peptides`. It is a 
##' `SingleCellExperiment` object containing the quantitative data for 13055 
##' peptides in 24 samples.
##' 
##' Sample annotation is stored in `colData(zhu2018MCP)`.
##' 
##' See `Details` for information about data collection.
##' 
##' @details 
##' 
##' **Acquisition protocole**
##' 
##' The data was acquired using the following setup. More information can be 
##' found in the source article (see `References`).
##' 
##' - **Cell isolation**: brain patches were collected using laser-captuer 
##'   microdissection (PALM MicroBeam) on flash frozen rat (*Rattus norvergicus*)
##'   brain tissues. Note that the samples were stained with H&E before 
##'   dissection for histological analysiis. DMSO is used as sample collection 
##'   solution
##' - **Sample preparation** performed using the nanoPOTs device. DMSO 
##'   evaporation + protein extraction (DMM + DTT) + alkylation (IAA) + Lys-C 
##'   digestion + trypsin digestion.
##' - **Separation**: nanoLC (Dionex UltiMate with an in-house packed 60cm x 
##'   30μm LC columns; 50nL/min)
##' - **Ionization**: ESI (2,000V)
##' - **Mass spectrometry**: Thermo Fisher Orbitrap Fusion Lumos Tribrid (MS1
##'   accumulation time = 246ms; MS1 resolution = 120,000; MS1 AGC = 3E6). The 
##'   MS/MS settings  depend on the sample size, excepted for the AGC = 1E5. 
##'   50µm (time = 502ms; resolution = 240,000), 100µm (time = 246ms; resolution 
##'   = 120,000), 200µm (time = 118ms; resolution = 60,000). 
##' - **Data analysis**: MaxQuant (v1.5.3.30) + Perseus (v1.5.6.0) + Origin Pro
##'   2017
##' 
##' **Data collection**
##' 
##' The PSM data were collected from the PRIDE repository (accession ID: 
##' PXD008844).  We downloaded the `MaxQuant_Peptides.txt` file containing the 
##' combined identification and quantification results. The sample annotation 
##' was infered from the names of columns holding the quantification data. The 
##' data were then converted to a `QFeatures` object using the `scp::readSCP` 
##' function. 
##' 
##' @source 
##' The PSM data can be downloaded from the PRIDE repository PXD008844.
##' 
##' @references 
##' Zhu, Ying, Maowei Dou, Paul D. Piehowski, Yiran Liang, Fangjun Wang, Rosalie 
##' K. Chu, William B. Chrisler, et al. 2018. “Spatially Resolved Proteome 
##' Mapping of Laser Capture Microdissected Tissue with Automated Sample 
##' Transfer to Nanodroplets.” Molecular & Cellular Proteomics: MCP 17 (9): 
##' 1864–74 ([link to article](http://dx.doi.org/10.1074/mcp.TIR118.000686)).
##' 
##' @docType data
##' 
##' @keywords datasets
##' 
##' 
"zhu2018MCP"

