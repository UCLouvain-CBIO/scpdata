##' Leduc, Andrew, Yanxin Xu, Gergana Shipkovenska, Zhixun Dou, and Nikolai
##' Slavov. 2025. “Limiting the Impact of Protein Leakage in Single-Cell
##' Proteomics.” Nature Communications 16 (1): 4169.
##'
##' @description
##'
##' Single-cell proteomics data from primary mouse tracheal epithelial cells
##' prepared fresh immediately after tissue dissociation. Each cell was stained
##' with Sytox Green before sorting, so the dataset carries a per-cell measure
##' of membrane permeabilisation alongside the protein quantifications. The
##' dataset contains quantitative information on 928 fresh single cells. Only
##' about 4% of these cells were permeabilised, which makes the dataset a
##' baseline for comparison with cells affected by protein leakage.
##'
##' Only the fresh cells are included. The frozen cells from the same study are
##' not part of this dataset because no identification file was available for
##' them.
##'
##' @format A [QFeatures] object with 34 sets, each being a [SingleCellExperiment]
##' object:
##'
##' - Sets 1-32: PSM data with identifications and quantitive data - results from
##'   a MaxQuant database search.
##' - Set 33 ("Peptides"): processed peptide data containing quantitative data coming
##'   from the PSM data. Taken from the source material.
##' - Set 44 ("Proteins"): processed protein data containing quantitative data coming
##'   from the peptide data.
##'
##' The `colData(leduc2025_pSCoPE())` contains cell type and batch annotations that
##' are common to all assays.
##'
##' @section Acquisition protocol:
##'
##' The data were acquired using the following setup. More information can be
##' found in the source article (see `References`).
##'
##' - **Tissue dissociation**: whole tracheas from 4-month-old male C57BL/6
##'   mice were dissociated for 30 minutes at 37 °C with a papain, collagenase I,
##'   hyaluronidase, DNase and dispase enzyme cocktail.
##' - **Cell isolation**: CellenOne sorting
##' - **Sample preparation**: nPOP sample preparation. Cell lysis (DMSO and heat)
##'   was followed by overnight trypsin digestion, then TMTpro 35-plex labelling.
##'   Labelled samples were pooled and dried for injection. Channel 126C was used
##'   as a carrier and 127N as a reference. Cells were labeled in sets of 29.
##' - **Separation**: online nLC (Vanquish Neo UHPLC) with a 25 cm x 75 µm
##'   IonOpticks column.
##' - **Mass spectrometry**: Orbitrap Exploris 480 with prioritised data
##'   acquisition (pSCoPE).
##' - **Data analysis**: MaxQuant 2.4.3.0 searched against the murine SwissProt
##'   database (downloaded 30 July 2018; 20,386 entries) at 1% FDR. Downstream
##'   analysis was performed in R using the QuantQC package.
##'
##' @section Data collection:
##'
##' Raw and searched mass spectrometry data are available from the MassIVE
##' repository (accession MSV000094790). Processed data are available from the
##' Slavov lab website (see `Source`).
##'
##' The necessary files for reproducing the QFeatures object are:
##'     - evidence.txt
##'     - celltype_annotation_fresh.csv
##'     - meta_fresh.csv
##'     - peptides.txt
##'     - proteinGroups.txt
##'
##' All data come from the author distribution.
##'
##' @source
##' The data were downloaded from the
##' [Slavov Lab](https://scp.slavovlab.net/Leduc_et_al_2024) website via a
##' shared Goodle Drive
##' [folder](https://drive.google.com/drive/folders/1IW3An_b7vKpfDbanjT2XqSt1mz0N6VUX).
##' The raw data, identification and quantification data can also be found in the
##' MassIVE repository `MSV000094790`:
##' ftp://massive-ftp.ucsd.edu/v07/MSV000094790/.
##'
##' The data used and generated as part of the scpdata repackaging is available
##' in the [Zenodo archive](https://zenodo.org/records/24079803).
##'
##' @references
##' Leduc, Andrew, Yanxin Xu, Gergana Shipkovenska, Zhixun Dou, and Nikolai
##' Slavov. 2025. “Limiting the Impact of Protein Leakage in Single-Cell
##' Proteomics.” Nature Communications 16 (1): 4169.
##' ([link to article](https://www.nature.com/articles/s41467-025-56736-7#Abs1)).
##'
##' @examples
##' \donttest{
##' leduc2025_pSCoPE()
##' }
##'
##' @keywords datasets
##'
##' @name leduc2025_pSCoPE
NULL
