
####---- SCPDATA PACKAGE MAN PAGE ----####


##' Single-Cell Proteomics Data Package
##' 
##' @description
##' 
##' The `scpdata` package distributes mass spectrometry-based 
##' single-cell proteomics datasets. The datasets were collected from 
##' published work and formatted to a standardized data framework. 
##' The `scp` frameworks stores the expression data for different MS 
##' levels (identified spectrum, peptide, or protein) in separate 
##' assays. Each assay is an object of class [SingleCellExperiment] 
##' that allows easy integration with state-of-the-art single-cell 
##' analysis tools. All assays are contained in a single object of 
##' class [QFeatures]. An overview of the data structure is shown 
##' provided in the `scp` package. 
##' 
##' The `scpdata()` function returns a summary table with all 
##' currently available datasets in the package. More information 
##' about the data content and the data collection can be found in the
##' corresponding manual pages. 
##' 
##' @return A `DataFrame` table containing a summary of the available 
##' datasets.
##' 
##' @seealso More information about the data manipulation can be found
##' in the `scp` package.
##' 
##' @examples 
##' ## List available datasets and their metadata 
##' scpdata()
##' 
##' ## Load data using the ExperimentHub interface
##' hub <- ExperimentHub()
##' 
##' \dontrun{
##' ## Download the data set of interest using ExperimentHub indexing
##' hub[["EH7711"]]
##' ## Download the same data set using the build-in function
##' leduc2022()
##' }
##' 
##' @author Christophe Vanderaa
##'
##' @aliases scpdata-package scpdata
##' 
##' @import ExperimentHub
##' @import QFeatures
##' @import SingleCellExperiment
##' @importFrom AnnotationHub query
##' @importFrom S4Vectors mcols
##' @export
##' 
scpdata <- function() {
    mcols(query(ExperimentHub(), "scpdata"))
}


#' Deprecated leduc2022 dataset
#'
#' The `leduc2022` dataset has been updated to include plexDIA and
#' pSCoPE data. The new datasets names are `leduc2022_pSCoPE`
#' (previously `leduc2022`) and `leduc2022_plexDIA` (new). See the
#' respective documentation pages for more information.
#'
#' @return The `leduc2022_pSCoPE` dataset.
#' 
#' @export
leduc2022 <- function() {
    warning(
        "Deprecated dataset! The leduc2022() dataset has been ", 
        "updated to include plexDIA and pSCoPE data. Available datasets ",
        "are now: leduc2022_pSCoPE (previously leduc2022) and ",
        "leduc2022_plexDIA (new). Returning leduc2022_pSCoPE()..."
    )
    leduc2022_pSCoPE()
}
