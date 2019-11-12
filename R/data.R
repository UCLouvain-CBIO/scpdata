####---- DESCRIPTION ----####

# Data utilities along with data documentation


#' List the available data sets in the scpdata package
#' 
#' The function lists the available datasets in the package. 
#' 
#' @import MSnbase
#' 
#' @details 
#' See the documentation of a particular data set for more information (eg 
#' \code{?specht2019}). 
#' 
#' @return 
#' An object of class "packageIQR" with available data sets.
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



####---- SPECHT ET AL. 2019 ----####


#' Quantifying the emergence of macrophage heterogeneity using the SCoPE2 
#' pipeline (Specht et al. 2019)
#'
#' Single cell proteomics data produced and published by Specht et al. from the 
#' Slavov Lab (see references). It contains quantitative information 356 cells. 
#' Cells can be either macrophages (n = 259) or monocytes (n = 97). 
#' 
#' @usage
#' data("specht2019_protein")
#' data("specht2019_peptide")
#' 
#' @format 
#' \code{specht2019} contains 2 sets of data:
#' \itemize{
#'   \item \code{specht2019_peptide}: an MSnSet with peptide expression levels 
#'   for 6787 peptides x 356 cells.
#'   \item \code{specht2019_protein}: an MSnSet with protein expression levels
#'   for 2316 proteins x 356 cells.
#' }
#' See Details for information about data collection.
#'
#' @details
#' 
#' Three data files are provided by Specht and colleagues are:
#' \itemize{
#'   \item \code{Peptides-raw.csv}: Peptides x single cells at 1% FDR. The 
#'   first 2 columns list the corresponding protein identifiers and peptide 
#'   sequences and each subsequent column corresponds to a single cell. Peptide 
#'   identification is based on spectra analyzed by MaxQuant and is enhanced by 
#'   using DART-ID to incorporate retention time information. 
#'   \item \code{Proteins-processed.csv}: Proteins x single cells at 1% FDR, 
#'   imputed and batch corrected.
#'   \item \code{Cells.csv}: Annotation x single cells. Each column corresponds 
#'   to a single cell and the rows include relevant metadata, such as, cell type 
#'   if known, measurements from the isolation of the cell, and derivative 
#'   quantities, i.e., rRI, CVs, reliability.
#'  }
#'  
#'  \strong{Peptide expression data: \code{specht2019_peptide}}
#'  
#'  The \code{Peptides-raw.csv} data have already been partially processed by 
#'  the authors. The MS files were analyzed with MaxQuant and DART-ID and the 
#'  output \code{evidence.txt} file was parsed into R. This output was further 
#'  processed by the authors as follows:
#'  \itemize{
#'    \item Keep only the single cell runs (experiments FP94 and FP97)
#'    \item Correct TMT reporter itensities (RI) for isotopic cross contamination.
#'    \item Remove out reverse hits and contaminants (identified by MaxQuant), 
#'    and contaminated spectra (\code{PIF > 0.8})
#'    \item Remove peptides with low identification score (\code{FDR >= 0.01} or 
#'    \code{PEP >= 0.02})
#'    \item Remove cells with less than 300 peptides
#'    \item Remove peptides that are more than 10\% the intensity of the carrier
#'    \item Divide peptide intensities in every channel by the reference channel
#'    \item Convert data to matrix format
#'    \item Zero or infinite intensities are replaced by \code{NA}'s
#'    \item Remove cells having \code{median CV > 0.43}, havuing 30th quantile of 
#'    the log10 transformed relative RIs smaller than -2.5, or having median of 
#'    the log10 transformed relative RI larger than -1.3
#'    \item Divide column (cells) with median intensity and divide rows with 
#'    mean intensity
#'    \item Remove rows (peptides) then columns (cells) that contain more than 
#'    99\% of missing data
#'    \item Log2 transform the data
#'  }
#'  
#'  The peptide data (\code{Peptides-raw.csv}) and the meta data 
#'  (\code{Cells.csv}) were combined into an \code{\link{MSnSet}} object.
#'  
#'  \strong{Protein expression data: \code{specht2019_protein}}
#'  
#'  On top of the steps described above, the peptide expression data were further
#'  processed by the authors to the protein data through the following steps: 
#'  \itemize{
#'    \item Aggregate the peptide to their corresponding protein using their 
#'    median value.
#'    \item Subtracting rows by the averages and subtract column by the column
#'    median
#'    \item Impute missing values using K-nearest neighbors (k = 3)
#'    \item Correct for batch effect using the \code{\link{ComBat}} algorithm.
#'  }
#'  This leads to the \code{Proteins-processed.csv} file. This protein 
#'  expression data are combined with the meta data (\code{Cells.csv}) into an 
#'  \code{\link{MSnSet}} object.#'  
#'  
#' @source 
#' The original data can be downloaded from the 
#' \href{https://scope2.slavovlab.net/docs/data}{Slavov Lab} website.
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
#' @aliases specht2019 specht2019_peptide
#'
"specht2019_protein"


####---- DOU ET AL. 2019 ----####


#' High-Throughput Single Cell Proteomics Enabled by Multiplex 
#' Isobaric Labeling in a Nanodroplet Sample Preparation Platform: 
#' HeLa digests (Dou et al. 2019)
#'  
#' @description 
#' Single-cell proteomics using nanoPOTS combined with TMT isobaric labeling. 
#' The data set is published as the supplementary data set 1 by Dou et al. in 
#' 2019 (see references). The expression matrix contains 1641 proteins x 20 
#' single-cell HeLa digests. The samples are not truly single cells but are 
#' commercial Hela digest diluted to single cell amounts (0.2ng). The boost 
#' wells contain the same digest but at hihgher dose (10 ng).
#' 
#' @usage 
#' data("dou2019_1_protein")
#' 
#' @format 
#' \code{dou2019_1} contains 1 set of data:
#' \itemize{
#'   \item \code{dou2019_1_protein}: an MSnSet with protein expression levels
#'   for 1641 proteins x 20 "cells".
#' }
#' See Details for information about data collection.
#' 
#' @details 
#' The data set was downloaded from the supplementary information section of the
#' publisher's website (see sources). 
#' 
#' \strong{Protein expression data: \code{dou2019_3_protein}}
#'
#' The spreadsheet is called \code{ac9b03349_si_003.xlsx} and contains 7 sheets 
#' from which we only took the sheet 6 (named \code{"5 - Run 1 and 2 raw data"}) 
#' containing the combined data of two MS runs. It corresponds to the assembly 
#' of MSGF+ identifications and MASIC reporter assemblies. The data is then isotope 
#' corrected and sum rolled-up to the protein level. Contaminant and reverse hit 
#' are removed from this table. This is all performed by the authors. We 
#' converted the data to an \code{\link{MSnSet}} object.
#' 
#' @source 
#' The original data can be downloaded from the 
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
#' @docType data
#' 
#' @keywords datasets
#' 
#' @aliases dou2019_1
#' 
"dou2019_1_protein"



#' High-Throughput Single Cell Proteomics Enabled by Multiplex 
#' Isobaric Labeling in a Nanodroplet Sample Preparation Platform: 
#' testing boosting ratios (Dou et al. 2019)
#'  
#' @description 
#' Single-cell proteomics using nanoPOTS combined with TMT isobaric labeling. 
#' The data set is published as the supplementary data set 2 by Dou et al. in 
#' 2019 (see references). The expression matrix contains 1436 proteins x 60 
#' single cells. The cell types are either "Raw" (macrophage cells), "C10" 
#' (epihelial cells), or "SVEC" (endothelial cells). Each cell was replicated 2 
#' or 3x. Each cell type was run using 3 levels of boosting: 0 ng (no boosting), 
#' 5 ng or 50 ng. When boosting was applied, 1 reference well and 1 boosting 
#' well were added, otherwise 1 empty well was added. Each boosting setting 
#' (no boosting, 5ng, 50ng) was run twice.
#' 
#' @usage 
#' data("dou2019_2_protein")
#' 
#' @format 
#' \code{dou2019_2} contains 1 set of data:
#' \itemize{
#'   \item \code{dou2019_2_protein}: an MSnSet with protein expression levels
#'   for 1436 proteins x 60 cells.
#' }
#' See Details for information about data collection.
#' 
#' @details 
#' The data set was downloaded from the supplementary information section of the
#' publisher's website (see sources). 
#' 
#' \strong{Protein expression data: \code{dou2019_2_protein}}
#'
#' The spreadsheet is called \code{ac9b03349_si_004.xlsx} and contains 7 sheets 
#' from which we took the 2nd, 4th and 6th sheets (named \code{"01 - No Boost 
#' raw data"}, \code{"03 - 5ng boost raw data"}, \code{"05 - 50ng boost raw 
#' data"}, respectively). It corresponds to the assembly of MSGF+ 
#' identifications and MASIC reporter assemblies. The data is then isotope 
#' corrected and sum rolled-up to the protein level. Contaminant and reverse hit 
#' are removed from this table. This is all performed by the authors. We matched 
#' the proteins between the three boosting settings and combined all data in a 
#' single table  The boosting quantities are kept as phenotype data. The data 
#' set is formated to an \code{\link{MSnSet}} object.
#' 
#' @source 
#' The original data can be downloaded from the 
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
#' @docType data
#' 
#' @keywords datasets
#' 
#' @aliases dou2019_2
#' 
"dou2019_2_protein"


#' High-Throughput Single Cell Proteomics Enabled by Multiplex 
#' Isobaric Labeling in a Nanodroplet Sample Preparation Platform: profiling of 
#' murine cell populations (Dou et al. 2019)
#'  
#' @description 
#' Single-cell proteomics using nanoPOTS combined with TMT isobaric labeling. 
#' The data set is published as the supplementary data set 3 by Dou et al. in 
#' 2019 (see references). The expression matrix contains 2331 proteins x 132 
#' wells. Among the 132 wells, 72 contained single cells, corresponding to 24 
#' C10 cells, 24 RAW cells, and 24 SVEC. The other wells are eithers boosting 
#' channels (12), empty channels (36) or reference channels (12). The 
#' different cell types where evenly distributed across 4 nanoPOTS chips. 
#' Samples were 11-plexed with TMT ions.
#' 
#' @usage 
#' data("dou2019_3_protein")
#' 
#' @format 
#' \code{dou2019_3} contains 1 set of data:
#' \itemize{
#'   \item \code{dou2019_3_protein}: an MSnSet with protein expression levels
#'   for 2331 proteins x 132 cells.
#' }
#' See Details for information about data collection.
#' 
#' @details 
#' The data set was downloaded from the supplementary information section of the
#' publisher's website (see sources). 
#' 
#' \strong{Protein expression data: \code{dou2019_3_protein}}
#'
#' The spreadsheet is called \code{ac9b03349_si_005.xlsx} and contains 7 sheets 
#' from which we took only the 2nd (named \code{"01 - Raw sc protein data"}). It 
#' corresponds to the assembly of MSGF+ identifications and MASIC reporter 
#' assemblies. The data is then isotope corrected and sum rolled-up to the 
#' protein level. Contaminant and reverse hit are removed from this table. This 
#' is performed by the authors. We converted the data to an \code{\link{MSnSet}} 
#' object.
#' 
#' @source 
#' The original data can be downloaded from the 
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
#' \code{\link{dou2019_1}}, \code{\link{dou2019_2}}
#' 
#' @docType data
#' 
#' @keywords datasets
#' 
#' @aliases dou2019_3
#' 
"dou2019_3_protein"


#' mPOP SCoPE-MS Master Mix 20180824 (Specht et al. 2018)
#'
#' @description 
#' Single cell proteomics data produced and published by Specht et al. from the 
#' Slavov Lab (see references). It contains quantitative information for 1824 
#' peptides x 220 channels. Channels are either carrier (50 cell equivalent of
#' Jurkat or U-937 digestion product), empty, or single cell equivalent of 
#' digestion product (Jurkat or U-937).
#' 
#' @usage 
#' data("specht2018_peptide")
#' 
#' @format 
#' \code{speccht2018} contains 1 set of data:
#' \itemize{
#'   \item \code{specht2018_peptide}: an MSnSet with peptide expression levels
#'   for 1838 peptides x 220 cells.
#' }
#' See Details for information about data collection.
#' 
#' @details
#' 
#' The peptide expression data was parsed from the \code{evidence.txt} file from
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
#'   \item All information is stored in an \code{\ling{MSnSet}} object
#' }
#'  We finally formated the data to an \code{\link{MSnSet}} object.
#'
#' @source 
#' The original data can be downloaded from the 
#' \href{http://slavovlab.net/mPOP/index.html}{Slavov Lab} website.
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
#' @aliases specht2018
#' 
"specht2018_peptide"