% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{zhu2018MCP}
\alias{zhu2018MCP}
\title{Zhu et al. 2018 (Mol. Cel. Prot.): rat brain laser dissections}
\format{
A \link{QFeatures} object with 1 assay called \code{peptides}. It is
a \link{SingleCellExperiment} object containing the quantitative data
for 13055 peptides in 24 samples.

Sample annotation is stored in \code{colData(zhu2018MCP())}.
}
\source{
The PSM data can be downloaded from the PRIDE repository
PXD008844. FTP link
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2018/07/PXD008844
}
\usage{
zhu2018MCP
}
\description{
Near single-cell proteomics data of laser captured
micro-dissection samples. The samples are 24 brain sections from
rat pups (day 17). The slices are 12 um thick squares of either
50, 100, or 200 um width. 5 samples were dissected from the corpus
callum (\code{CC}), 4 samples were dissected from the corpus collosum
(\code{CP}), 13 samples were extracted from the cerebral cortex
(\code{CTX}), and 2 samples are labeled as (\code{Mix}).
}
\section{Acquisition protocol}{


The data were acquired using the following setup. More information
can be found in the original article (see \code{References}).
\itemize{
\item \strong{Cell isolation}: brain patches were collected using
laser-capture microdissection (PALM MicroBeam) on flash frozen
rat (\emph{Rattus norvergicus}) brain tissues. Note that the samples
were stained with H&E before dissection for histological
analysis. DMSO is used as sample collection solution
\item \strong{Sample preparation} performed using the nanoPOTs device: DMSO
evaporation + protein extraction (DMM + DTT) + alkylation (IAA)
\itemize{
\item Lys-C digestion + trypsin digestion.
}
\item \strong{Separation}: nanoLC (Dionex UltiMate with an in-house packed
60cm x 30um LC columns; 50nL/min)
\item \strong{Ionization}: ESI (2,000V)
\item \strong{Mass spectrometry}: Thermo Fisher Orbitrap Fusion Lumos
Tribrid (MS1 accumulation time = 246ms; MS1 resolution =
120,000; MS1 AGC = 3E6). The MS/MS settings  depend on the
sample size, excepted for the AGC = 1E5. 50um (time = 502ms;
resolution = 240,000), 100um (time = 246ms; resolution =
120,000), 200um (time = 118ms; resolution = 60,000).
\item \strong{Data analysis}: MaxQuant (v1.5.3.30) + Perseus (v1.5.6.0) +
Origin Pro 2017
}
}

\section{Data collection}{


The PSM data were collected from the PRIDE repository (accession
ID: PXD008844).  We downloaded the \code{MaxQuant_Peptides.txt} file
containing the combined identification and quantification results.
The sample annotation was inferred from the names of columns
holding the quantification data. The data were then converted to a
\link{QFeatures} object using the \link[scp:readSCP]{scp::readSCP} function.
}

\examples{
zhu2018MCP()

}
\references{
Zhu, Ying, Maowei Dou, Paul D. Piehowski, Yiran Liang, Fangjun
Wang, Rosalie K. Chu, William B. Chrisler, et al. 2018. “Spatially
Resolved Proteome Mapping of Laser Capture Microdissected Tissue
with Automated Sample Transfer to Nanodroplets.” Molecular &
Cellular Proteomics: MCP 17 (9): 1864–74
(\href{http://dx.doi.org/10.1074/mcp.TIR118.000686}{link to article}).
}
\keyword{datasets}
