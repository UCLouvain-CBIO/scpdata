####---- DESCRIPTION ----####

# Data utilities along with data documentation


##' List the available data sets in the `scpdata` package
##' 
##' The function lists the available datasets in the package. 
##' 
##' @import QFeatures
##' @import SingleCellExperiment
##' 
##' @details 
##' See the documentation of a particular data set for more information (eg 
##' `?specht2019v2`). 
##' 
##' @return 
##' An object of class `"packageIQR"` containing all available data sets.
##' 
##' @examples
##' scpdata()
##' 
##' @export
##' 
scpdata <- function(){
  
  out <- data(package = "scpdata")
  return(out)
}

####---- SPECHT ET AL. 2019 ----####


##' FACS + SCoPE2: macrophages vs monocytes (Specht et al. 2019 - version 2)
##'
##' Single cell proteomics data acquired by the Slavov Lab. This is the version 
##' 2 of the data released in December 2019. It contains quantitative 
##' information of macrophages and monocytes at PSM, peptide and protein level. 
##' 
##' @format A `Features` object with 179 assays, each assay being a 
##' `SingleCellExperiment` object: 
##' \itemize{
##'   \item Assay 1-63: PSM data for SCoPE2 sets acquired with a TMT-11 
##'   multiplexing protocole, hence those assays contain 11 columns. Columns 
##'   hold quantitative information from single-cell channels, carrier channels, 
##'   reference channels, empty (blank) channels and unused channels.
##'   \item Assay 64-177: PSM data for SCoPE2 sets acquired with a TMT-16
##'   multiplexing protocole, hence those assays contain 16 columns. Columns 
##'   hold quantitative information from single-cell channels, carrier channels, 
##'   reference channels, empty (blank) channels and unused channels.
##'   \item Assay 178: peptide data containing quantitative data for 9208 
##'   peptides and 1018 single-cells. Cell type annotation and batch annotation
##'   are stored in `colData(specht2019v2[[178]]`.
##'   \item Assay 179: peptide data containing quantitative data for 2772 
##'   proteins and 1018 single-cells. Cell type annotation and batch annotation
##'   are stored in `colData(specht2019v2[[179]]`.
##' }
##' The `colData(specht2019v2)` contains cell type annotation and batch 
##' annotation for all assays.
##' 
##' See `Details`` for information about data collection.
##'
##' @details The PSM data were collected from a shared Goolge Drive folder that 
##' is accessible from the SlavovLab website. The folder contains the following 
##' files of interest: 
##' \itemize{
##'   \item `ev_updated.txt`: the MaxQuant/DART-ID output file
##'   \item `annotation_fp60-97.csv`: sample annotation
##'   \item `batch_fp60-97.csv`: batch annotation
##' }
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
##' \href{https://scope2.slavovlab.net/docs/data}{Slavov Lab} website via a
##' shared Google Drive 
##' \href{https://drive.google.com/drive/folders/1VzBfmNxziRYqayx3SP-cOe2gu129Obgx}{folder}. 
##' The raw data and quantification data can also be found in the massIVE 
##' repository 
##' \href{ftp://massive.ucsd.edu/MSV000083945}{MSV000083945}.
##' 
##' @references Specht, Harrison, Edward Emmott, Toni Koller, and Nikolai Slavov. 
##' 2019. “High-Throughput Single-Cell Proteomics Quantifies the Emergence of 
##' Macrophage Heterogeneity.” bioRxiv (\href{https://doi.org/10.1101/665307}{DOI}).
##' 
##' @docType data
##'
##' @keywords datasets
##' 
"specht2019v2"