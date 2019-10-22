####---- DESCRIPTION ----####

# Data utilities along with data documentation


#' Available data sets
#' 
#' Function listing the available data sets in the \code{scpdata} package
#'
#' @export
#'
scpdata <- function(){
  data(package = "scpdata")
}

####---- SPECHT ET AL. 2019 ----####

#' 356 cells x 6787 peptides data containing macrophages and monocytes (Specht et al. 2019)
#'
#' Single cell proteomics data produced and published by Specht et al. from the 
#' Slavov Lab (see references). It contains quantitative information for 6787 
#' peptides in 356 cells. Cells can be either macrophages (n = 259) or monocytes
#' (n = 97). The data is formated to an MSnSet object with no further 
#' modification of the original data. 
#' 
#' \strong{Title}: High-throughput single-cell proteomics quantifies the emergence of macrophage heterogeneity
#' 
#' \strong{Abstract}: The fate and physiology of individual cells are controlled by networks of 
#' proteins. Yet, our ability to quantitatively analyze protein networks in 
#' single cells has remained limited. To overcome this barrier, we developed 
#' SCoPE2. It integrates concepts from Single-Cell ProtEomics by Mass 
#' Spectrometry (SCoPE-MS) with automated and miniaturized sample preparation, 
#' substantially lowering cost and hands-on time. SCoPE2 uses data-driven 
#' analytics to optimize instrument parameters for sampling more ion copies per 
#' protein, thus supporting quantification with improved count statistics. 
#' These advances enabled us to analyze the emergence of cellular heterogeneity 
#' as homogeneous monocytes differentiated into macrophage-like cells in the 
#' absence of polarizing cytokines. We used SCoPE2 to quantify over 2,000 
#' proteins in 356 single monocytes and macrophages in about 85 hours of 
#' instrument time, and the quantified proteins allowed us to discern single 
#' cells by cell type. Furthermore, the data uncovered a continuous gradient of 
#' proteome states for the macrophage-like cells, suggesting that macrophage 
#' heterogeneity may emerge even in the absence of polarizing cytokines. Our 
#' methodology lays the foundation for quantitative analysis of protein networks 
#' at single-cell resolution.
#'
#' @docType data
#'
#' @usage data("specht2019")
#'
#' @format An object of class \code{"\link{MSnSet}"}
#'
#' @keywords datasets
#'
#' @references Specht, Harrison, Edward Emmott, Toni Koller, and Nikolai Slavov. 
#' 2019. “High-Throughput Single-Cell Proteomics Quantifies the Emergence of 
#' Macrophage Heterogeneity.” bioRxiv. https://doi.org/10.1101/665307.
#' (\href{https://www.biorxiv.org/content/10.1101/665307v2}{bioRxiv})
#'
#' @source The original data can be downloaded from the 
#' \href{https://scope2.slavovlab.net/docs/data}{Slavov Lab} website.
#'
"specht2019"
