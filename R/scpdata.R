
####---- SCPDATA PACKAGE MAN PAGE ----####


##' Single Cell Proteomics Data Package
##' 
##' Data package distributing mass spectrometry-based single-cell proteomics 
##' datasets from published work. The datasets are curated and formated to 
##' a standardized data framework. This frameworks stores the different assays
##' (different batches for PSM data, peptide data and protein data) in 
##' `SingleCellExperiment` objects. The assays are stored in a single 
##' `QFeatures` object. 
##' 
##' @seealso More information about the data framework can be found in the 
##' [scp] package
##' 
##' 
##' @examples 
##' library(ExperimentHub)
##' 
##' ## Load data using the ExperimentHub interface
##' hub <- ExperimentHub()
##' ## x <- query(hub, c("scpdata", "specht2019v2"))
##' 
##' \dontrun{
##' ## download resource
##' dataset <- x[[1]]
##' }
##'
##' ## Load data by name
##' x <- specht2019v2(metadata = FALSE)
##' 
##' 
##' @aliases scpdata-package scpdata
"scpdata-package"


####---- SCPDATA UTILITIES ----####


#' @import ExperimentHub
#' 
#' Get `scpdata` metadata
#' 
#' The function returns the metadata table of the datasets that are available in
#' `scpdata`.
#'
#' @return A `DataFrame` table containing the dataset metadata.
#' 
#' @references 
#' See in the respective data sets' manual pages for references to publications.
#' 
#' @export
#'
#' @examples
#' meta <- scpdata()
#' ## Number of available datasets
#' nrow(meta)
#' ## List datasets 
#' meta$Title
#' 
#' 
scpdata <- function() {
  mcols(query(ExperimentHub, "scpdata"))
}
