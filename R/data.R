####---- DESCRIPTION ----####

# Data utilities along with data documentation


#' Available data sets
#' 
#' Function listing the available data sets in the \code{scpdata} package
#'
#' @export
#'
scpdata <- function(){
  out <- data(package = "scpdata")
  return(out)
}


####---- SPECHT ET AL. 2019 ----####


#'  Quantifying the emergence of macrophage heterogeneity using the SCoPE2 
#'  pipeline (Specht et al. 2019)
#'
#' Single cell proteomics data produced and published by Specht et al. from the 
#' Slavov Lab (see references). It contains quantitative information for 6787 
#' peptides x 356 cells. Cells can be either macrophages (n = 259) or monocytes
#' (n = 97). The data are formated to an MSnSet object. 
#' 
#' @details
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
#' Macrophage Heterogeneity.” bioRxiv (\href{https://doi.org/10.1101/665307}{DOI}).
#'
#' @source The original data can be downloaded from the 
  #' \href{https://scope2.slavovlab.net/docs/data}{Slavov Lab} website.
#'
"specht2019"


####---- DOU ET AL. 2019 ----####


#' Deep proteome coverage for single cell analysis using nanoPOTS combined with
#' TMT isobaric labeling method (Dou et al. 2019)
#' 
#' @description 
#' Single cell proteomics data produced and published by Dou et al. (see 
#' references). The article comes with 3 supplementary datasets: \code{dou2019_1},
#' \code{dou2019_2}, and \code{dou2019_3} (see Details for more information).
#' All data are formated to an MSnSet object. 
#' 
#' @aliases dou2019_1 dou2019_2 dou2019_3
#' 
#' @details  
#' \strong{Title}: High-Throughput Single Cell Proteomics Enabled by Multiplex 
#' Isobaric Labeling in a Nanodroplet Sample Preparation Platform.
#' 
#' \strong{Abstract}:  Effective extension of mass spectrometry-based proteomics 
#' to single cells remains challenging. Herein we combined microfluidic 
#' nanodroplet technology with tandem mass tag (TMT) isobaric labeling to 
#' significantly improve analysis throughput and proteome coverage for single 
#' mammalian cells. Isobaric labeling facilitated multiplex analysis of single 
#' cell-sized protein quantities to a depth of ∼1 600 proteins with a median CV 
#' of 10.9% and correlation coefficient of 0.98. To demonstrate in-depth high 
#' throughput single cell analysis, the platform was applied to measure protein 
#' expression in 72 single cells from three murine cell populations (epithelial, 
#' immune, and endothelial cells) in <2 days instrument time with over 2 300 
#' proteins identified. Principal component analysis grouped the single cells 
#' into three distinct populations based on protein expression with each 
#' population characterized by well-known cell-type specific markers. Our 
#' platform enables high throughput and unbiased characterization of single cell 
#' heterogeneity at the proteome level.
#' 
#' \strong{\code{dou2019_1}}: Supplementary data set 1, raw data for HeLa 
#' digest. The expression matrix contains 1641 proteins x 20 single-cell digests. Note that the samples are not truly 
#' single cells but are commercial Hela digest diluted to single cell amounts
#' (0.2ng). The boost wells contain the same digest but at hihgher dose (10 ng).
#' 
#' \strong{\code{dou2019_2}}: Supplementary data set 2, raw data for testing 
#' boosting ratios. The expression matrix contains 1436 proteins x 60 single 
#' cells. The cell type are either \code{"Raw"} (macrophage cells), \code{"C10"}
#' (epihelial cells), or \code{"SVEC"} (endothelial cells). Each cell was 
#' replicated 2 or 3x. When boosting (5ng or 50ng) was applied, 1 reference 
#' well and 1 boosting well were added, otherwise (no boosting) 1 empty well 
#' was added. Each boosting setting (no boosting, 5ng, 50ng) was run twice.
#' 
#' \strong{\code{dou2019_3}}: Supplementary data set 3, raw data for isobaric 
#' labelling-based single cell quantification and bulk-scale label free 
#' quantification. The expression matrix contains 2331 proteins x 132 wells. Among 
#' the 132 well, 72 contained single cells, corresponding to 24 C10 cells, 24 RAW cells, 
#' and 24 SVEC. The other wells are eithers boosting channels (12), empty channels
#' (36) or reference channels (12). Note that the different cell types where
#' evenly distributed across 4 nanoPOTS chips. Samples were 11-plexed with
#' TMT ions. 
#' 
#' \strong{\code{dou2019_3}}: Supplementary data set 3, raw data for isobaric 
#' labelling-based single cell quantification and bulk-scale label free 
#' quantification 
#' 
#' @docType data
#' 
#' @usage 
#' dou2019()
#' data("dou2019_1")
#' data("dou2019_2")
#' data("dou2019_3")
#' 
#' @format  
#' \code{dou2019()} loads all three data sets at once. Each data set
#' is an object of class \code{"\link{MSnSet}"}.
#' 
#' @keywords datasets
#' 
#' @references Dou, Maowei, Geremy Clair, Chia-Feng Tsai, Kerui Xu, William B. 
#' Chrisler, Ryan L. Sontag, Rui Zhao, et al. 2019. “High-Throughput Single Cell 
#' Proteomics Enabled by Multiplex Isobaric Labeling in a Nanodroplet Sample 
#' Preparation Platform.” Analytical Chemistry, September 
#' (\href{https://doi.org/10.1021/acs.analchem.9b03349}{DOI}).
#' 
#' @source The original data can be downloaded from the \href{https://pubs.acs.org/doi/10.1021/acs.analchem.9b03349}{ACS Publications}
#' website (Supplementary information section).
#' 
#' @export
dou2019 <- function(){
  data("dou2019_1")
  data("dou2019_2")
  data("dou2019_3")
}


#'  mPOP SCoPE-MS Master Mix 20180824 (Specht et al. 2018)
#'
#' Single cell proteomics data produced and published by Specht et al. from the 
#' Slavov Lab (see references). It contains quantitative information for 1824 
#' peptides x 220 channels. Channels are either carrier (50 cell equivalent of
#' Jurkat or U-937 digestion product), empty, or single cell equivalent of 
#' digestion product (Jurkat or U-937). The data are formated to an MSnSet object. 
#' 
#' @details
#' 
#' \strong{Title}: Automated sample preparation for high-throughputsingle-cell 
#' proteomics
#' 
#' \strong{Abstract}: A major limitation to applying quantitative LC-MS/MS 
#' proteomics to small samples, such as single cells, are the losses incured 
#' during sample cleanup. To relieve this limitation, we developed a Minimal 
#' ProteOmic sample Preparation (mPOP) method for culture-grown mammalian 
#' cells. mPOP obviates cleanup and thus eliminates cleanup-related losses while 
#' expediting sample preparation and simplifying its automation.  Bulk SILAC 
#' samples processed by mPOP or by conventional urea-based methods indicated 
#' that mPOP results in complete cell lysis and accurate relative 
#' quantification. We integrated mPOP lysis with the Single CellProtEomics by 
#' Mass Spectrometry (SCoPE-MS) sample preparation, and benchmarked the 
#' quantification of such samples on a Q-exactive instrument. The results 
#' demonstrate low noise and high technical reproducibility. Then, we FACS 
#' sorted single U-937, HEK-293, andmouse ES cells into 96-well plates and 
#' analyzed them by automated mPOP and SCoPE-MS.The quantified proteins enabled 
#' separating the single cells by cell-type and cell-division-cycle phase.
#'
#' @docType data
#'
#' @usage data("specht2018")
#'
#' @format An object of class \code{"\link{MSnSet}"}
#'
#' @keywords datasets
#'
#' @references 
#' Specht, Harrison, Guillaume Harmange, David H. Perlman, Edward Emmott, 
#' Zachary Niziolek, Bogdan Budnik, and Nikolai Slavov. 2018. “Automated Sample 
#' Preparation for High-Throughput Single-Cell Proteomics.” bioRxiv 
#' (\href{https://doi.org/10.1101/399774}{DOI}).
#'
#' @source The original data can be downloaded from the 
#' \href{http://slavovlab.net/mPOP/index.html}{Slavov Lab} website.
#'
"specht2018"