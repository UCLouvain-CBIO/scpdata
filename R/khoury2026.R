##' Khoury et al, 2026 (bioRxiv): Single-cell proteomics reveals
##' cell-type-specific functional coordination in PBMCs
##'
##' @description
##'
##' Single-cell samples were prepared using the nPOP glass-slide sample
##' preparation method and multiplexed with 2-plexDIA (mTRAQ). Proteomics data
##' were acquired with slice-PASEF on timsTOF Ultra, timsTOF Ultra 2 and
##' timsUltra AIP mass spectrometers. The dataset contains quantitative
##' information on 1,651 single peripheral blood mononuclear cells
##' (monocytes, NK cells, CD8 T cells, CD4 T cells and B cells) from two healthy
##' donors. The data are available at the PSM, peptide and protein levels.
##'
##' The paper investigates within-cell-type protein covariation (protein
##' complex coordination, protein stability, functional rewiring between cell
##' types and cytokine covariation) using matching single-cell RNA sequencing
##' (Smart-seq3xpress) data from the same donors.
##'
##' @format A [QFeatures] object with 1103 assays, each assay being a
##' [SingleCellExperiment] object:
##'
##' - Assay 1-1097: DIA-NN main output report table. Each run acquires 1 or 2
##' single cells. It contains the results of the spectrum identification and
##' quantification.
##' - `Peptide_Normalized_Centered`: peptide data containing quantitative data
##' for 23,572 peptides and 1,651 single-cells. Values are column normalized
##' and centered, so they are relative.
##' - `Peptide_Column_Normalized`: peptide data containing quantitative data
##' for 23,572 peptides and 1,651 single-cells. Values are column normalized and
##' log transformed, so they are absolute values.
##' - `Protein_Imputed_Column_Normalized`: protein data containing quantitative
##' data for 3,732 proteins and 1,651 single-cells with k-nearest neighbors (KNN)
##' imputation and comBat batch-correction.
##' - `Protein_Imputed_Normalized_Centered`: protein data containing
##' quantitative data for 3,732 proteins and 1,651 single-cells with k-nearest
##' neighbors (KNN) imputation and comBat batch-correction.
##' - `Protein_Normalized_Centered`: protein data containing quantitative data
##' for 3,732 proteins and 1,651 single-cells.
##' - `Protein_Column_Normalized`: protein data containing quantitative data
##' for 3,732 proteins and 1,651 single-cells.
##'
##' The `colData(khoury2026())` contains cell type and batch annotations that
##' are common to all assays.
##'
##' @section Acquisition protocol:
##'
##' The data were acquired using the following setup. More information
##' can be found in the source article (see `References`).
##'
##' - **Cell isolation**: CellenONE cell sorting.
##' - **Sample preparation** performed using the nPOP sample preparation method.
##' cell lysis (DMSO) + trypsin digestion + mTRAQ labelling (d0/d4) + pooling
##' with an mTRAQ d8-labelled isotopologous carrier (5 cell-equivalents).
##' - **Separation**: online nLC (NanoElute2, UltiMate 3000 or Vanquish Neo
##' HPLC) with a 25cm x 75um IonOpticks Aurora column; 30 min active gradient
##' (4-40% buffer B) at 200 nL/min.
##' - **Mass spectrometry**: timsTOF Ultra, timsTOF Ultra 2, and timsUltra AIP
##' mass spectrometers.
##' - **Data analysis**: DIA-NN (v1.8.2 beta 17).
##'
##' @section Data collection:
##'
##' The data were collected from a Zenodo folder that is accessible from
##' the SlavovLab website (see `Source` section).
##'
##' The folder ('Batch_IDs_and_abundances/') contains the following
##' files of interest:
##' - `IdentificationsList_n2p1.csv`: PSM level processed quantitative data.
##' - `IdentificationsList_n2p2.csv`: PSM level processed quantitative data.
##' - `IdentificationsList_n3p1.csv`: PSM level processed quantitative data.
##' - `IdentificationsList_n6p1.csv`: PSM level processed quantitative data.
##' - `IdentificationsList_n6p2.csv`: PSM level processed quantitative data.
##'
##' Following files of interest also available in Zenodo:
##' - `Peptide_x_Cell_Matrix_Normalized_Centered.csv`: peptide level processed
##' quantitative data.
##' - `Peptide_x_Cell_Matrix_Column_Normalized.csv`: peptide level processed
##' quantitative data.
##' - `Protein_x_Cell_Matrix_Imputed_Column_Normalized.csv`: protein level
##' processed quantitative data.
##' - `Protein_x_Cell_Matrix_Imputed_Normalized_Centered.csv`: protein level
##' processed quantitative data.
##' - `Protein_x_Cell_Matrix_Normalized_Centered.csv`: protein level processed
##' quantitative data.
##' - `Protein_x_Cell_Matrix_Column_Normalized.csv`: protein level processed
##' quantitative data.
##' - `PBMC_covariation_protein_cell_metadata.csv`: Sample metadata.
##' - `PepProtMap_integrated.csv`: Peptide to Protein mapping.
##'
##' We formatted the PSM level quantification tables so that columns match
##' with those of the annotation and peptide/protein level data. Each PSM level
##' data turned into [SingleCellExperiment] objects including row data
##' information. All the PSM level assays used to generate [QFeatures] object.
##'
##' The peptide and protein assays were formatted as [SingleCellExperiment]
##' objects and the sample metadata were matched to the column names. The
##' objects are then added to the [QFeatures] object and the rows of the PSM
##' data are linked to the rows of the corresponding peptide data based on the
##' sequence information through an `AssayLink` object. Peptide assays and
##' protein assays linked in the same way through Gene Symbols as provided by
##' the authors.
##'
##' @source
##' The data were downloaded from the
##' [Slavov Lab](https://scp.slavovlab.net/Khoury_et_al_2026) website via a
##' shared Zenodo folder
##' [folder](https://zenodo.org/records/22649483).
##'
##' The data used and generated as part of the scpdata repackaging is available
##' in the [Zenodo archive](https://doi.org/10.5281/zenodo.24007418).
##'
##' @references
##' Khoury L, et al., Single-cell proteomics reveals cell-type-specific
##' functional coordination in PBMCs, bioRxiv, doi: 10.64898/2026.09.11.751068
##' ([link to article](https://www.biorxiv.org/content/10.64898/2026.09.11.751068v1)).
##'
##' @examples
##' \donttest{
##' khoury2026()
##' }
##'
##' @keywords datasets
##'
##' @name koury2026
NULL
