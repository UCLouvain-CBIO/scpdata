####---- DESCRIPTION ----####

# Data utilities along with data documentation


#' List the available data sets in the scpdata package
#' 
#' The function lists the available datasets in the package. 
#' 
#' @import MSnbase
#' @import SingleCellExperiment
#' 
#' @details 
#' See the documentation of a particular data set for more information (eg 
#' \code{?specht2019}). 
#' 
#' @return 
#' An object of class "packageIQR" containing all available data sets.
#' 
#' @examples
#' scpdata()
#' 
#' @export
#' 
scpdata <- function(){
  out <- data(package = "scpdata")
  return(out)
}

####---- ZHU ET AL. 2018, MOL. AND CELL. PROTEOMICS ----####

#' laser dissection + nanoPOTS:  (Zhu et al. 2018, Mol Cell Proteomics)
#' 
#' Near single-cell proteomics data produced and published by Zhu et al. 
#' (see references). It contains quantitative information for 24 mouse brain 
#' sections. The samples were brains from rat pups (day 17). The slices are 12
#' µm thick squares of either 50, 100, or 200 µm width. 5 samples were dissected 
#' from the corpus callum (\code{"CC"}), 4 samples were dissected from the 
#' corpus collosum (\code{"CP"}), 13 samples were extracted from the cerebral 
#' cortex (\code{"CTX"}), and 2 samples are labeled as (\code{"Mix"}).
#' 
#' @usage
#' data("zhu2018MCP_peptide")
#' 
#' @format 
#' \itemize{
#'   \item \code{zhu2018MCP_peptide}: a \code{\link{SingleCellExperiment}} 
#'   with peptide expression levels for 10074 peptides x 24 samples
#' }
#' See Details for information about data collection.
#' 
#' @details
#' 
#' Peptide data were collected from the PRIDE repository (accession ID: 
#' PXD008844), sample annotations were infered from column names.
#' 
#' \strong{Peptide expression data: \code{zhu2018MCP_peptide}}
#' 
#' We processed the MaxQuant output \code{MaxQuant_Peptides.txt} as follows:
#' \itemize{
#'   \item Extract the phenotype data (sample type) from the intensity column 
#'   names.
#'   \item Remove contaminants and reverse hits 
#'   \item Remove peptides with a posterior error probability < 0.01
#'   \item Replace 0's by NA's
#'   \item Remove rows (peptides) that contain only NA's
#' }
#' We finally formated the data to an \code{\link{SingleCellExperiment}} 
#' object. The samples in this data set are coming from different brain regions
#' (CTX, CP, or CC) and dissected as 50, 100 or 200 µm wide squares.
#' 
#' @source 
#' The raw data, the idenfication and quantification results can be found in the
#' PRIDE repository 
#' \href{ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2018/07/PXD008844}{PXD008844}
#' 
#' @references 
#' Zhu, Ying, Maowei Dou, Paul D. Piehowski, Yiran Liang, Fangjun Wang, Rosalie 
#' K. Chu, William B. Chrisler, et al. 2018. “Spatially Resolved Proteome 
#' Mapping of Laser Capture Microdissected Tissue with Automated Sample
#' Transfer to Nanodroplets.” Molecular & Cellular Proteomics: MCP 17 (9): 
#' 1864–74. \href{http://dx.doi.org/10.1074/mcp.TIR118.000686}{DOI}
#' 
#' @docType data
#'
#' @keywords datasets
#' 
#' @aliases zhu2018MCP_peptide
#' 
"zhu2018MCP_peptide"

####---- ZHU ET AL. 2018, NATURE COMM ----####

#' laser dissection + nanoPOTS: T1D islets (Zhu et al. 2018, Nat Comm)
#' 
#' Near single-cell proteomics data produced and published by Zhu et al. 
#' (see references). It contains quantitative information for 18 human pancreas 
#' samples. The samples are pancreatic islets laser dissected from pancreatic 
#' tissue slices from pancreata recovered from organ donors through the 
#' JDRFNetwork for Pancreatic Organ Donors with Diabetes (nPOD) program. 9 
#' samples out 18 are control islets, the 9 others are from type 1 diabetes 
#' (T1D) patients. 
#' 
#' @usage
#' data("zhu2018NC_islets_peptide")
#' 
#' @format 
#' \itemize{
#'   \item \code{zhu2018NC_islets_peptide}: a \code{\link{SingleCellExperiment}} 
#'   with peptide expression levels for 18611 peptides x 18 samples
#' }
#' See Details for information about data collection.
#' 
#' @details
#' 
#' Peptide data were collected from the PRIDE repository (accession ID: 
#' PXD006847), sample annotations were infered from column names.
#' 
#' \strong{Peptide expression data: \code{zhu2018NC_islets_peptide}}
#' 
#' We processed the MaxQuant output \code{Islet_t1d_ct_peptides.txt} as 
#' follows:
#' \itemize{
#'   \item Extract the phenotype data (sample type) from the intensity column 
#'   names.
#'   \item Remove contaminants and reverse hits 
#'   \item Remove peptides with a posterior error probability < 0.01
#'   \item Replace 0's by NA's
#' }
#' We finally formated the data to an \code{\link{SingleCellExperiment}} 
#' object. The samples in this data set are either a cell lysate 
#' (\code{"lysate"}), HeLA cells (\code{"cells"}), blank samples 
#' (\code{"blank"}), or cancer cell lines (\code{"MCF7"} and \code{"THP1"}).
#' 
#' @source 
#' The raw data, the idenfication and quantification results can be found in the
#' PRIDE repository 
#' \href{ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2018/01/PXD006847}{PXD006847}
#' 
#' @references 
#' Zhu, Ying, Paul D. Piehowski, Rui Zhao, Jing Chen, Yufeng Shen, Ronald J. 
#' Moore, Anil K. Shukla, et al. 2018. “Nanodroplet Processing Platform for Deep 
#' and Quantitative Proteome Profiling of 10-100 Mammalian Cells.” Nature 
#' Communications 9 (1): 882. \href{http://dx.doi.org/10.1038/s41467-018-03367-w}{DOI}
#' 
#' @docType data
#'
#' @keywords datasets
#' 
#' @aliases zhu2018NC_islets_peptide
#' 
"zhu2018NC_islets_peptide"


#' nanoPOTS: HeLa dilutions (Zhu et al. 2018, Nat Comm)
#' 
#' Near single-cell proteomics data produced and published by Zhu et al. 
#' (see references). It contains quantitative information for 21 samples, 
#' some being triplicates of HeLa dilution (~10, ~40 and ~140 cells). Remaining
#' samples are blanks (0 cells), lysate or cancer cell lines (THP-1 or MCF-7).
#' 
#' @usage
#' data("zhu2018NC_hela_peptide")
#' 
#' @format 
#' \itemize{
#'   \item \code{zhu2018NC_hela_peptide}: a \code{\link{SingleCellExperiment}} 
#'   with peptide expression levels for 30345 peptides x 21 samples
#' }
#' See Details for information about data collection.
#' 
#' @details
#' 
#' Peptide data were collected from the PRIDE repository (accession ID: 
#' PXD006847), sample annotations were infered from column names.
#' 
#' \strong{Peptide expression data: \code{zhu2018NC_hela_peptide}}
#' 
#' We processed the MaxQuant output \code{CulturedCells_peptides.txt} as 
#' follows:
#' \itemize{
#'   \item Extract the phenotype data (sample type, and number of cells when 
#'   single cells) from the intensity column names.
#'   \item Remove contaminants and reverse hits 
#'   \item Remove peptides with a posterior error probability < 0.01
#'   \item Replace 0's by NA's
#' }
#' We finally formated the data to an \code{\link{SingleCellExperiment}} 
#' object. The samples in this data set are either a cell lysate 
#' (\code{"lysate"}), HeLA cells (\code{"cells"}), blank samples 
#' (\code{"blank"}), or cancer cell lines (\code{"MCF7"} and \code{"THP1"}).
#' 
#' @source 
#' The raw data, the idenfication and quantification results can be found in the
#' PRIDE repository 
#' \href{ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2018/01/PXD006847}{PXD006847}
#' 
#' @references 
#' Zhu, Ying, Paul D. Piehowski, Rui Zhao, Jing Chen, Yufeng Shen, Ronald J. 
#' Moore, Anil K. Shukla, et al. 2018. “Nanodroplet Processing Platform for Deep 
#' and Quantitative Proteome Profiling of 10-100 Mammalian Cells.” Nature 
#' Communications 9 (1): 882. \href{http://dx.doi.org/10.1038/s41467-018-03367-w}{DOI}
#' 
#' @docType data
#'
#' @keywords datasets
#' 
#' @aliases zhu2018NC_hela_peptide
#' 
"zhu2018NC_hela_peptide"


####---- SPECHT ET AL. 2019 ----####


#' FACS + SCoPE2: macrophages vs monocytes (Specht et al. 2019)
#'
#' Single cell proteomics data produced and published by Specht et al. from the 
#' Slavov Lab (see references). It contains quantitative information 356 cells. 
#' Cells can be either macrophages (n = 259) or monocytes (n = 97). 
#' 
#' @usage
#' data("specht2019v1_protein")
#' data("specht2019v1_peptide")
#' 
#' @format 
#' \itemize{
#'   \item \code{specht2019v1_peptide}: a \code{\link{SingleCellExperiment}} 
#'   with peptide expression levels for 7515 peptides x 682 samples
#'   \item \code{specht2019v1_protein}: a \code{\link{SingleCellExperiment}} 
#'   with protein expression levels for 2316 proteins x 356 samples (=cells)
#' }
#' See Details for information about data collection.
#'
#' @details
#' 
#' Peptide data were collected from the massIVE database (accession ID: 
#' MSV000083945), annotations were received during a personnal discussion but 
#' are also contained in version 2 of the data (see 
#' \code{\link{specht2019v2_peptide}}): 
#' \itemize{
#'   \item \code{ev_updated.txt}: the MaxQuant/DART-ID output file
#'   \item \code{annotation_fp60-97.csv}: sample phenotype annotation
#'   \item \code{batch_fp60-97.csv}: batch annotation
#' }
#'  
#' Protein data were collected from the author's website (see sources).
#' \itemize{
#'   \item \code{Proteins-processed.csv}: Proteins x single cells at 1% FDR, 
#'   imputed and batch corrected.
#'   \item \code{Cells.csv}: Annotation x single cells. Each column corresponds 
#'   to a single cell and the rows include relevant metadata, such as, cell type 
#'   if known, measurements from the isolation of the cell, and derivative 
#'   quantities, i.e., rRI, CVs, reliability.
#' }
#'  
#' \strong{Peptide expression data: \code{specht2019v1_peptide}}
#'  
#' We processed the MaxQuant output \code{ev_updated.txt} as follows:
#' \itemize{
#'   \item Keep only column of interest, that is the peptide and protein 
#'   information, the quantification information (TMT intensities), the 
#'   identification information (PEP, q-values, PIF,...) and the contamination
#'   information
#'   \item Keep only the single cell runs (experiments FP94 and FP97)
#'   \item Remove out reverse hits and contaminants (identified by MaxQuant), 
#'   and contaminated spectra (\code{PIF > 0.8})
#'   \item Remove peptides with low identification score (\code{FDR >= 0.01} or 
#'   \code{PEP >= 0.02})
#'   \item Remove cells with less than 300 peptides
#'   \item Some peptides were identified twice within the same experiment set 
#'   because of different charge states. The peptide identification with lowest 
#'   PEP is kept and the other(s) discarded.
#'   \item Intensity data is formated to a feature (peptide) x sample 
#'   (combination of experiment run and TMT channel) matrix.
#'   \item Zero intensities are replaced by NA's
#'   \item Peptide information is gathered in a feature data frame, only 
#'   peptide information common to all runs is kept (eg sequence, mass, 
#'   protein,...).
#'   \item Sample information is gathered in a phenotype data frame. The sample 
#'   information is extracted from the \code{annotation_fp60-97.csv} and 
#'   \code{batch_fp60-97.csv} files.
#' }
#' 
#' We finally formated the data to an \code{\link{SingleCellExperiment}} 
#' object. The samples in this data set are either a carrier well 
#' (\code{"carrier_mix"}), a normalization well (\code{"norm"}), an unused
#' well (\code{"unused"}), an empty well (\code{"sc_0"}), a macrophage cell
#' (\code{"sc_m0"}), or a monocyte cell (\code{"sc_u"}).
#' 
#' \strong{Protein expression data: \code{specht2019v1_protein}}
#'  
#' On top of the steps described above, the peptide expression data were further
#' processed by the authors to the protein data through the following steps: 
#' \itemize{
#'   \item Remove peptides that are more than 10\% the intensity of the carrier
#'   \item Divide peptide intensities in every channel by the reference channel
#'   \item Zero or infinite intensities are replaced by NA's
#'   \item Remove cells having median CV > 0.43, having 30th quantile of the 
#'   log10 transformed relative RIs smaller than -2.5, or having median of the 
#'   log10 transformed relative RI larger than -1.3
#'   \item Divide column (cells) with median intensity and divide rows with mean 
#'   intensity
#'   \item Remove rows (peptides) then columns (cells) that contain more than 
#'   99\% of missing data
#'   \item Log2 transform the data
#'   \item Aggregate the peptide to their corresponding protein using their 
#'   median value.
#'   \item Normalize rows by subtracting the row average and normalize columns 
#'   by subtracting the column median
#'   \item Impute missing values using K-nearest neighbors (k = 3)
#'   \item Correct for batch effect using the \code{\link{ComBat}} algorithm.
#' }
#' 
#' This leads to the \code{Proteins-processed.csv} file. The protein 
#' expression data are combined with the meta data (\code{Cells.csv}) into an 
#' \code{\link{SingleCellExperiment}} object. The data set contains 
#' macrophages (\code{"sc_m0"}) and monocytes (\code{"sc_u"}). 
#'  
#' @source 
#' The data can be downloaded from the 
#' \href{https://scope2.slavovlab.net/docs/data}{Slavov Lab} website. The raw 
#' data and quantification data can also be found in the massIVE repository 
#' \href{ftp://massive.ucsd.edu/MSV000083945}{MSV000083945}
#' 
#' @references 
#' Specht, Harrison, Edward Emmott, Toni Koller, and Nikolai Slavov. 
#' 2019. “High-Throughput Single-Cell Proteomics Quantifies the Emergence of 
#' Macrophage Heterogeneity.” bioRxiv (\href{https://doi.org/10.1101/665307}{DOI}).
#' 
#' @docType data
#'
#' @keywords datasets
#' 
#' @seealso 
#' \code{\link{specht2019v2_peptide}}, 
#' \code{\link{specht2019v2_protein}}
#' 
#' @aliases specht2019v1_protein specht2019v2_peptide
#'
"specht2019v1_peptide"


####---- DOU ET AL. 2019 ----####


#' FACS + nanoPOTS + TMT: HeLa digests (Dou et al. 2019)
#'  
#' @description 
#' Single-cell proteomics using nanoPOTS combined with TMT isobaric labeling. 
#' The samples are not truly single cells but are commercial Hela digest diluted 
#' to single cell amounts (0.2ng). The boosting wells contain the same digest 
#' but at higher dose (10 ng).
#' 
#' @usage 
#' data("dou2019_hela_protein")
#' data("dou2019_hela_peptide")
#' 
#' @format 
#' \itemize{
#'   \item \code{dou2019_hela_protein}: a \code{\link{SingleCellExperiment}}
#'   object with protein expression levels for 24,797 peptides x 20 "cells".
#'   \item \code{dou2019_hela_protein}: a \code{\link{SingleCellExperiment}}
#'   object with protein expression levels for 1,641 proteins x 20 "cells".
#'   
#' }
#' See Details for information about data collection.
#' 
#' @details 
#' 
#' \strong{Peptide expression data: \code{dou2019_hela_peptide}}
#' 
#' The peptide data was constructed as follows. For every MS run:
#' \itemize{
#'   \item Load the MSGF+ identification files 
#'   (\code{"Hela_run_*_msgfplus.mzid"}) and remove decoy and contaminant PSMs
#'   \item Load the MASIC quantification files 
#'   (\code{"Hela_run_*_ReporterIons.txt"})
#'   \item Combine the identification and the quantification data by matching 
#'   the \code{"scan.number.s"} and \code{"scanNumber"} fields, respectively
#'   \item Format the data as an \code{\link{MSnSet}} object keeping 
#'   \item Perform TMT-10 istope correction 
#'   \item Sum-roll the PSM intensities to peptide intensities
#'   \item Keep only useful feature data (peptide and protein sequence, start 
#'   and end position of peptide location, protein length, protein description)
#'   \item Extract the sample information (TMT ion, run number and experiment 
#'   label) from the files names. Furthermore, sample type (lysate, carrier or
#'   blank) and sample amount (in ng) were added from the article (see 
#'   Reference).
#' }
#' Next, combine all runs in a single MSnSet (matching on peptide name), replace
#' NA's by 0's and convert to an \code{\link{SingleCellExperiment}} object.
#' 
#' \strong{Protein expression data: \code{dou2019_hela_protein}}
#'
#' The data set is published as the supplementary data set 1 by Dou et al. in 
#' 2019 (see references). The expression matrix contains 1641 proteins x 20 
#' HeLa digests. The spreadsheet is called \code{ac9b03349_si_003.xlsx} and 
#' contains 7 sheets from which we only took the sheet 6 (named \code{"5 - Run 
#' 1 and 2 raw data"}) containing the combined data of two MS runs. It 
#' corresponds to the assembly of MSGF+ identifications and MASIC reporter 
#' assemblies. The data is then isotope corrected and sum rolled-up to the 
#' protein level. Contaminant and reverse hit are removed from this table. This 
#' is all performed by the authors. We converted the data to a 
#' \code{\link{SingleCellExperiment}} object.
#' 
#' @source 
#' The peptide and raw data can be downloaded from the massIVE repository 
#' \href{ftp://massive.ucsd.edu/MSV000084110/}{MSV000084110}.
#' 
#' The protein data can be downloaded from the 
#' \href{https://pubs.acs.org/doi/10.1021/acs.analchem.9b03349}{ACS Publications}
#' website (Supplementary information section).
#' 
#' @references 
#' Dou, Maowei, Geremy Clair, Chia-Feng Tsai, Kerui Xu, William B. 
#' Chrisler, Ryan L. Sontag, Rui Zhao, et al. 2019. “High-Throughput Single Cell 
#' Proteomics Enabled by Multiplex Isobaric Labeling in a Nanodroplet Sample 
#' Preparation Platform.” Analytical Chemistry, September 
#' (\href{https://doi.org/10.1021/acs.analchem.9b03349}{DOI}).
#' 
#' @seealso 
#' \code{\link{dou2019_mouse_protein}}, \code{\link{dou2019_boosting_protein}}
#' 
#' @docType data
#' 
#' @keywords datasets
#' 
#' @aliases dou2019_hela_peptide dou2019_hela_protein
#' 
"dou2019_hela_protein"



#' FACS + nanoPOTS + TMT: testing boosting ratios (Dou et al. 2019)
#'  
#' @description 
#' Single-cell proteomics using nanoPOTS combined with TMT isobaric labeling. 
#' The cell types are either "Raw" (macrophage cells), "C10" 
#' (epihelial cells), or "SVEC" (endothelial cells). Each cell is replicated 2 
#' or 3x. Each cell type was run using 3 levels of boosting: 0 ng (no boosting), 
#' 5 ng or 50 ng. When boosting was applied, 1 reference well and 1 boosting 
#' well were added, otherwise 1 empty well was added. Each boosting setting 
#' (no boosting, 5ng, 50ng) was run twice.
#' 
#' @usage
#' data("dou2019_boosting_protein")
#' data("dou2019_boosting_peptide")
#' 
#' @format 
#' \itemize{
#'   \item \code{dou2019_boosting_peptide}: a \code{\link{SingleCellExperiment}}
#'    with protein expression levels for 51,567 peptides x 60 sample.
#'   \item \code{dou2019_boosting_protein}: a \code{\link{SingleCellExperiment}}
#'    with protein expression levels for 1,436 proteins x 60 samples.
#' }
#' See Details for information about data collection.
#' 
#' @details 
#' 
#' \strong{Peptide expression data: \code{dou2019_boosting_peptide}}
#' 
#' The peptide data was constructed as follows. For every MS run:
#' \itemize{
#'   \item Load the MSGF+ identification files 
#'   (\code{"Boosting_*_msgfplus.mzid"}) and remove decoy and 
#'   contaminant PSMs
#'   \item Load the MASIC quantification files 
#'   (\code{"Boosting_*_ReporterIons.txt"})
#'   \item Combine the identification and the quantification data by matching 
#'   the \code{"scan.number.s"} and \code{"scanNumber"} fields, respectively
#'   \item Format the data as an \code{\link{MSnSet}} object keeping 
#'   \item Perform TMT-10 istope correction 
#'   \item Sum-roll the PSM intensities to peptide intensities
#'   \item Keep only useful feature data (peptide and protein sequence, start 
#'   and end position of peptide location, protein length, protein description)
#'   \item Extract the sample information (TMT ion, run, boosting amount and 
#'   experiment label) from the files names. Furthermore, sample type was 
#'   obtained from Table S2 (see reference).
#' }
#' Next, combine all runs in a single MSnSet (matching on peptide name), replace
#' NA's by 0's and convert to an \code{\link{SingleCellExperiment}} object.
#' 
#' 
#' \strong{Protein expression data: \code{dou2019_boosting_protein}}
#'
#' The data set is published as the supplementary data set 2 by Dou et al. in 
#' 2019 (see references). The expression matrix contains 1436 proteins x 60 
#' single cells. The spreadsheet is called \code{ac9b03349_si_004.xlsx} and contains 7 sheets 
#' from which we took the 2nd, 4th and 6th sheets (named \code{"01 - No Boost 
#' raw data"}, \code{"03 - 5ng boost raw data"}, \code{"05 - 50ng boost raw 
#' data"}, respectively). It corresponds to the assembly of MSGF+ 
#' identifications and MASIC reporter assemblies. The data is then isotope 
#' corrected and sum rolled-up to the protein level. Contaminant and reverse hit 
#' are removed from this table. This is all performed by the authors. We matched 
#' the proteins between the three boosting settings and combined all data in a 
#' single table. The boosting quantities are kept as phenotype data. The data 
#' set is formated to an \code{\link{SingleCellExperiment}} object.
#' 
#' @source 
#' The peptide and raw data can be downloaded from the massIVE repository 
#' \href{ftp://massive.ucsd.edu/MSV000084110/}{MSV000084110}.
#' 
#' The protein data can be downloaded from the 
#' \href{https://pubs.acs.org/doi/10.1021/acs.analchem.9b03349}{ACS Publications}
#' website (Supplementary information section).
#' 
#' @references 
#' Dou, Maowei, Geremy Clair, Chia-Feng Tsai, Kerui Xu, William B. 
#' Chrisler, Ryan L. Sontag, Rui Zhao, et al. 2019. “High-Throughput Single Cell 
#' Proteomics Enabled by Multiplex Isobaric Labeling in a Nanodroplet Sample 
#' Preparation Platform.” Analytical Chemistry, September 
#' (\href{https://doi.org/10.1021/acs.analchem.9b03349}{DOI}).
#' 
#' @seealso 
#' \code{\link{dou2019_hela_protein}}, \code{\link{dou2019_mouse_protein}}
#' 
#' @docType data
#' 
#' @keywords datasets
#' 
#' @aliases dou2019_boosting_protein dou2019_boosting_peptide
#' 
"dou2019_boosting_protein"


#' FACS + nanoPOTS + TMT: profiling of murine cell populations (Dou et al. 2019)
#'  
#' @description 
#' Single-cell proteomics using nanoPOTS combined with TMT isobaric labeling. 
#' The cell types are either "Raw" (macrophage cells), "C10" (epihelial cells), 
#' or "SVEC" (endothelial cells). Out of the 132 wells, 72 contain single cells, 
#' corresponding to 24 C10 cells, 24 RAW cells, and 24 SVEC. The other wells are 
#' either boosting channels (12), empty channels (36) or reference channels 
#' (12). Boosting and reference channels are balanced (1:1:1) mixes of C10, 
#' SVEC, and RAW samples at 5 ng and 0.2 ng, respectively. The different cell 
#' types where evenly distributed across 4 nanoPOTS chips. Samples were 
#' 11-plexed with TMT labeling.
#' 
#' @usage 
#' data("dou2019_mouse_peptide")
#' data("dou2019_mouse_protein")
#'  
#' @format 
#' \itemize{
#'   \item \code{dou2019_mouse_peptide}: a \code{\link{SingleCellExperiment}}
#'   with protein expression levels for 74,374 proteins x 132 samples.
#'   \item \code{dou2019_mouse_protein}: a \code{\link{SingleCellExperiment}}
#'   with protein expression levels for 2331 proteins x 132 samples.
#' }
#' See Details for information about data collection.
#' 
#' @details 
#' 
#' \strong{Peptide expression data: \code{dou2019_mouse_peptide}}
#' 
#' The peptide data was constructed as follows. For every MS run:
#' \itemize{
#'   \item Load the MSGF+ identification files 
#'   (\code{"Single_Cell_Chip_*_msgfplus.mzid"}) and remove decoy and 
#'   contaminant PSMs
#'   \item Load the MASIC quantification files 
#'   (\code{"Single_Cell_Chip_*_ReporterIons.txt"})
#'   \item Combine the identification and the quantification data by matching 
#'   the \code{"scan.number.s"} and \code{"scanNumber"} fields, respectively
#'   \item Format the data as an \code{\link{MSnSet}} object keeping 
#'   \item Perform TMT-11 istope correction 
#'   \item Sum-roll the PSM intensities to peptide intensities
#'   \item Keep only useful feature data (peptide and protein sequence, start 
#'   and end position of peptide location, protein length, protein description)
#'   \item Extract the sample information (TMT ion, chip number and experiment 
#'   label) from the files names. Furthermore, sample type was obtained from 
#'   Table S1 (see reference).
#' }
#' Next, combine all runs in a single MSnSet (matching on peptide name), replace
#' NA's by 0's and convert to an \code{\link{SingleCellExperiment}} object.
#' 
#' \strong{Protein expression data: \code{dou2019_mouse_protein}}
#'
#' The data set is published as the supplementary data set 3 by Dou et al. in 
#' 2019 (see references). The expression matrix contains 2331 proteins x 132 
#' wells. The spreadsheet is called \code{ac9b03349_si_005.xlsx} and contains 7 
#' sheets from which we took only the 2nd (named \code{"01 - Raw sc protein 
#' data"}). It corresponds to the assembly of MSGF+ identifications and MASIC 
#' reporter assemblies. The data is then isotope corrected and sum rolled-up to 
#' the protein level. Contaminant and reverse hit are removed from this table. 
#' This is performed by the authors. We converted the data to a
#' \code{\link{SingleCellExperiment}} object.
#' 
#' @source 
#' The peptide and raw data can be downloaded from the massIVE repository 
#' \href{ftp://massive.ucsd.edu/MSV000084110/}{MSV000084110}.
#' 
#' The protein data can be downloaded from the 
#' \href{https://pubs.acs.org/doi/10.1021/acs.analchem.9b03349}{ACS Publications}
#' website (Supplementary information section).
#' 
#' @references 
#' Dou, Maowei, Geremy Clair, Chia-Feng Tsai, Kerui Xu, William B. 
#' Chrisler, Ryan L. Sontag, Rui Zhao, et al. 2019. “High-Throughput Single Cell 
#' Proteomics Enabled by Multiplex Isobaric Labeling in a Nanodroplet Sample 
#' Preparation Platform.” Analytical Chemistry, September 
#' (\href{https://doi.org/10.1021/acs.analchem.9b03349}{DOI}).
#' 
#' @seealso 
#' \code{\link{dou2019_hela_protein}}, \code{\link{dou2019_boosting_protein}}
#' 
#' @docType data
#' 
#' @keywords datasets
#' 
#' @aliases dou2019_mouse_protein dou2019_mouse_peptide
#' 
"dou2019_mouse_protein"


#' SCoPE-MS + mPOP lysis upgrade: Master Mix 20180824 (Specht et al. 2018)
#'
#' @description 
#' Single cell proteomics data produced and published by Specht et al. from the 
#' Slavov Lab (see references). Channels are either carrier (50 cell equivalent 
#' of Jurkat or U-937 digestion product), empty, or single cell equivalent of 
#' digestion product (Jurkat or U-937).
#' 
#' @usage 
#' data("specht2018_peptide")
#' 
#' @format 
#' \itemize{
#'   \item \code{specht2018_peptide}: a \code{\link{SingleCellExperiment}} 
#'   with peptide expression levels for 1838 peptides x 220 cells.
#' }
#' See Details for information about data collection.
#' 
#' @details
#' 
#' \strong{Peptide expression data: \code{specht2018_peptide}}
#' 
#' It contains quantitative information for 1824 peptides x 220 channels. The 
#' peptide expression data was parsed from the \code{evidence.txt} file from
#' the MaxQuant output. We processed the file as follows: 
#' \itemize{
#'   \item Contaminants and reverse hits are removed
#'   \item Peptides with PIF (parent ion fraction) smaller than 0.8 and PEP 
#'   (posterior error probability) greater then 0.02 are removed
#'   \item Experiment sets with less than 300 identified peptides are discarded. 
#'   Note that all experiment sets contained more than 300 identified peptides.
#'   \item Some peptides were identified twice within the same experiment set 
#'   because of different charge states. The peptide identification with lowest 
#'   PEP is kept and the other(s) discarded. 
#'   \item Intensity data is formated to a feature (peptide) x sample 
#'   (combination of experiment run and TMT channel) matrix. 
#'   \item Peptide information is gathered in a feature data frame. PEP, PIF, 
#'   Score, and retention time are aggregated by taking the median. When the 
#'   mass is differing among the same peptide sequence, only the mass for the 
#'   unmodified peptide is kept.
#'   \item Sample information is gathered in a phenotype data frame.
#' }
#'  We finally formated the data to a \code{\link{SingleCellExperiment}} object.
#'
#' \strong{Protein expression data: \code{specht2018_protein}}
#' 
#' TODO
#' 
#' @source 
#' The peptide data can be downloaded from the 
#' \href{http://slavovlab.net/mPOP/index.html}{Slavov Lab} website. The raw 
#' data and quantification data can also be found in the massIVE repository 
#' \href{ftp://massive.ucsd.edu/MSV000082841}{MSV000082841}
#'
#' @references 
#' Specht, Harrison, Guillaume Harmange, David H. Perlman, Edward Emmott, 
#' Zachary Niziolek, Bogdan Budnik, and Nikolai Slavov. 2018. “Automated Sample 
#' Preparation for High-Throughput Single-Cell Proteomics.” bioRxiv 
#' (\href{https://doi.org/10.1101/399774}{DOI}).
#' 
#' 
#' @docType data
#'
#' @keywords datasets
#' 
#' @aliases specht2018_peptide
#' 
"specht2018_peptide"