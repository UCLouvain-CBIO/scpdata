##' Leduc, Andrew, Yanxin Xu, Gergana Shipkovenska, Zhixun Dou, and Nikolai
##' Slavov. 2025. “Limiting the Impact of Protein Leakage in Single-Cell
##' Proteomics.” Nature Communications 16 (1): 4169.
##'
##' @description
##'
##' Samples were prepared using the nPOP sample preparation method for
##' multiplexed single cell proteomics1.Samples were run on an Exploris 480 mass
##' spectrometer with a Van- quish Neo liquid chromatography and autosampler. A
##' 25 cm 75 μmID IonOptics column was used for the chromatography column. The
##' gradient ramped from 8% to 40% buffer B (80% Acetonitrile 20% 0.1% formic
##' acid) over the course of 28 min with a 4 min wash at 90% buffer B at the
##' end. Mass spectrometry data acquisition of single cells was performed using
##' MaxQuant Live for Prioritized data acquisition with 60k MS1 and MS2
##' resolution, 118 ms maximum injection time and 1e6 maximum AGC.
##'
##' I stopped here ############################################
##'
##' Single-cell samples were prepared using the nPOP glass-slide sample
##' preparation method. Proteomics data were acquired using the timsControl 3.1
##' on a timsTOF SCP mass spectrometer. The dataset contains quantitative
##' information on 4,184 single cell from murine trachea. The data are available
##' at the peptide, and protein levels. The paper investigates the protein
##' regulation by translation and protein clearance using a matching single
##' cell RNA Sequencing data.
##'
##' @format A [QFeatures] object with 11 assays, each assay being a
##' [SingleCellExperiment] object:
##'
##' - Assay 1-4: unprocessed peptide data containing quantitative data coming
##'   from four replicates.
##' - Assay 5-8: processed peptide data containing quantitative data coming
##'   from four replicates.
##' - `proteins_absolute`: protein data containing quantitative data for 2,339
##'   proteins and 4,184 single-cells.
##' - `proteins_relative`: protein data containing quantitative data for 2,331
##'   proteins and 4,184 single-cells. Relative abundance values against the
##'   protein mean across all samples.
##' - `proteins_imputed`: protein data containing quantitative data for 667
##'   proteins and 4,184 single-cells with k-nearest neighbors (KNN) imputation
##'   and comBat batch-correction.
##'
##' The `colData(leduc2025())` contains cell type and batch annotations that
##' are common to all assays.
##'
##' @section Acquisition protocol:
##'
##' The data were acquired using the following setup. More information
##' can be found in the source article (see `References`).
##'
##' - **Cell isolation**: CellenONE cell sorting.
##' - **Sample preparation** performed using the nPOP sample preperation method.
##'   cell lysis (DMSO) + trypsin digestion + mTRAQ reagents.
##' - **Separation**: online nLC (Vanquish Neo UHPLC with a 25cm x 75um
##'   IonOpticks Aurora Series UHPLC column (AUR2-25075C18A).
##' - **Mass spectrometry**: timsTOF SCP mass spectrometer (MS1 scan range =
##'   100-1700 m/z; MS2 scan range = 300-1000 m/z).
##' - **Data analysis**: DIA-NN(v1.9.0).
##'
##' @section Data collection:
##'
##' The Peptide data were collected from a Zenodo folder that is accessible from
##' the SlavovLab website (see `Source` section).
##'
##' The folder ('02_raw_reptide_X_singleCell/') contains the following
##' files of interest:
##' - `r1_peptide.csv`: peptide level unprocessed quantitative data for
##' replicate 1.
##' - `r2_peptide.csv`: peptide level unprocessed quantitative data for
##' replicate 2.
##' - `r3_peptide.csv`: peptide level unprocessed quantitative data for
##' replicate 3.
##' - `r4_peptide.csv`: peptide level unprocessed quantitative data for
##' replicate 4.
##'
##' ('03_QuantQC_objects/') contains the following files of interest:
##' - `r1_5day_male.RData`: peptide level processed quantitative data for
##' replicate 1 and row data information for peptide to protein mapping.
##' - `r2_5day_female.RData`: peptide level processed quantitative data for
##' replicate 2 and row data information for peptide to protein mapping.
##' - `r3_10day_male.RData`: peptide level processed quantitative data for
##' replicate 3 and row data information for peptide to protein mapping.
##' - `r4_10day_female.RData`: peptide level processed quantitative data for
##' replicate 4 and row data information for peptide to protein mapping.
##'
##' ('04_Gene_X_SingleCell_and_annotations/') contains the following files of
##' interest:
##' `sc_protein_absolute.csv`: protein level quantitative data with absolute
##' protein concentrations.
##' `sc_protein_relative.csv`: protein level quantitative data with relative
##' protein abundance.
##' `sc_protein_annotations.csv`: sample and batch annotations.
##'
##' We also obtained imputed and batch corrected protein level data. To do that
##' we run the script `get_protein_Data.R` based on the original script
##' `02_cell_type_assign.R` in the github repository of the article.
##'
##' We formatted the peptide level quantification tables so that columns match
##' with those of the annotation and protein level data. Each peptide level
##' data for replicates turned into [SingleCellExperiment] objects including
##' row data information that is extracted from shared QuantQC objects. All
##' the peptide level assays used to generate [QFeatures] object.
##'
##' Protein level data were taken from '04_Gene_X_SingleCell_and_annotations/'
##' folder. Imputed and batch effect corrected protein level data is generated
##' using `get_protein_Data.R` script. The protein assays were formatted as
##' [SingleCellExperiment] objects and the sample metadata were matched to the
##' column names. The objects are then added to the [QFeatures] object
##' and the rows of the peptide data are linked to the rows of the corresponding
##' protein data based on the protein sequence information through an
##' `AssayLink` object.
##'
##' @source
##' The data were downloaded from the
##' [Slavov Lab](https://scp.slavovlab.net/Leduc_et_al_2025) website via a
##' shared Zenodo folder
##' [folder](https://zenodo.org/records/14902834).
##' The raw data and the quantification data can also be found in the
##' MassIVE repository `MSV000098940`:
##' ftp://massive-ftp.ucsd.edu/v10/MSV000098940/.
##'
##' The data used and generated as part of the scpdata repackaging is available
##' in the [Zenodo archive](https://doi.org/10.5281/zenodo.21904283).
##'
##' @references
##' Leduc A, et al., Principles of protein abundance regulation across single
##' cells in a mammalian tissue, bioRxiv, doi: 10.1101/2025.09.17.676955
##' ([link to article](https://www.biorxiv.org/content/10.1101/2025.09.17.676955v2)).
##'
##' @examples
##' \donttest{
##' leduc2025()
##' }
##'
##' @keywords datasets
##'
##' @name leduc2025
NULL
