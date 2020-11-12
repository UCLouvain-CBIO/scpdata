
####---- SCPDATA PACKAGE MAN PAGE ----####


##' Single-Cell Proteomics Data Package
##' 
##' @description
##' 
##' The `scpdata` package distributes mass spectrometry-based single-cell proteomics 
##' datasets. The datasets were collected from published work and formatted to 
##' a standardized data framework. The `scp` frameworks stores the expression 
##' data for different MS levels (identified spectrum, peptide, or protein) in 
##' separate assays. Each assay is an object of class `SingleCellExperiment` 
##' that allows easy integration with state-of-the-art single-cell analysis 
##' tools. All assays are contained in a single object of class `QFeatures`.
##' An overview of the data structure is shown below:
##' 
##' \figure{SCP_framework.pdf}
##' 
##' The `scpdata()` function returns a summary table with all currently 
##' available datasets in the package. More information about the data content 
##' and the data collection can be found in the corresponding manual pages. 
##' 
##' @return A `DataFrame` table containing a summary of the available datasets.
##' 
##' @seealso More information about the data manipulation can be found in the 
##' [scp] package.
##' 
##' @examples 
##' ## List available datasets and their metadata 
##' scpdata()
##' 
##' ## Load data using the ExperimentHub interface
##' hub <- ExperimentHub()
##' 
##' \dontrun{
##' ## download the data set of interest using EH indexing
##' hub[["EH3899"]]
##' ## download the same data set using scpdata function
##' `Specht et al. 2019 - SCoPE2 (biorRxiv): macrophages vs monocytes (version 2)`()
##' }
##'
##' @author Christophe Vanderaa
##'
##' @aliases scpdata-package scpdata
##' 
##' @import ExperimentHub
##' @import scp
##' 
scpdata <- function() {
  mcols(query(ExperimentHub(), "scpdata"))
}
