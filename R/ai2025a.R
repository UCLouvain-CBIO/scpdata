##' Ai et al., Cardiomyocytes Single Cell Proteomics, MCP (2025)
##'
##' @description
##'
##' The manuscript analysed cardiomyocytes derived from iPSC cells and adult
##' hearts. This annotated `QFeatures` objects contains the complete data, as
##' analysed with the `scp` package. The re-analysis was performed as part of
##' the EuBIC 2025 developer meeting hackathon dedicated to single-cell
##' proteomics data (https://lgatto.github.io/2025-EuBIC-SCP-hackathon/).
##'
##' The data contains a total of 304 assay. The 299 first assays contain the
##' individual single cardiomyocytes acquistions. The following assays contain
##' the joined precursor data, the log-transformed precuror data, the peptide
##' data, the protein data, and the scplainer-modelled precursor data (see the
##' link above for details).
##'
##' The `QFeatures` object was generated from the DIA-NN report file, available
##' from the `MsDataHub` package (see `?MsDataHub::Ai2025_aCMs_report.tsv`). The
##' report file is downloaded from the MassIVE MSV000094438
##' (doi:10.25345/C5T727S7Q) repository.
##'
##' @references
##'
##' - Ai, Lizhuo, Aleksandra Binek, Vladimir Zhemkov, Jae Hyung Cho, Ali
##' Haghani, Simion Kreimer, Edo Israely, et al. 2025. “Single Cell Proteomics
##' Reveals Specific Cellular Subtypes in Cardiomyocytes Derived from Human
##' iPSCs and Adult Hearts.” Molecular & Cellular Proteomics: MCP, no. 100910
##' (January): 100910.
##'
##' @keywords dataset
##'
##' @examples
##'
##' \donttest{
##' ai2025a()
##' }
"ai2025a"