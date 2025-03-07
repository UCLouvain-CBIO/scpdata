##' Hu et al, 2023 (The Journal of Physical Chemistry B): Correlated protein
##' modules
##'
##' @description
##'
##' They demonstrate the correlations between the levels of pairs of proteins
##' in single-cell proteomics (SCP) at steady state. In measuring pairwise
##' correlations among 1000 proteins in a population of K562 cells and oocytes,
##' they observed many correlated protein modules (CPMs) that are functionally
##' involved in certain biological functions. Certain CPMs are specific to a
##' particular cell type, some common to different cell types. Additionally,
##' compared to single-cell transcriptomics and bulk proteomics,
##' protein correlations are functionally and experimentally more significant
##' in SCP than those corresponding mRNAs.
##'
##' @format Two [SingleCellExperiment] objects:
##'
##' - `proteins_K562`: protein data containing quantitative data for 1249
##'   proteins and 69 single-cells with zero imputation.
##' - `proteins_oocyte`: protein data containing quantitative data for 3422
##'   proteins and 137 single-cells with zero imputation.
##'
##' The `colData(hu2023_oocyte())` contains cell type annotation.
##' The `colData(hu2023_K562())` contains cell type annotation.
##'
##' @section Acquisition protocol:
##'
##' The data were acquired using the following setup. More information
##' can be found in the source article (see `References`).
##'
##' - **Cell isolation**: K562 cells were re-suspended and washed in cold PBS.
##'   Single cells/10 cells were sorted into 96-well plates using a FACSAria
##'   instrument. Oocyte-cumulus complexes from C57/6J mice were collected
##'   after PMSG and HCG injections, with hyaluronidase used to remove cumulus
##'   cells. All samples stored at -80 degrees Celsius.
##' - **Sample preparation** Cells were digested with trypsin at 37 degrees
##'   Celsius for 3 hours. For label-free proteomics, digestion was terminated
##'   by adding 0.43% TFA and 1% ACN in water, followed by drying in a
##'   concentrator. Peptides were resuspended in 0.1% TFA and 1% ACN, and
##'   then transferred to sample tubes for LC-MS/MS analysis.
##' - **Separation**: 4 microliters of peptide digests were injected into a
##'   high-performance chromatography column (IonOpticks) and separated at a
##'   flow rate of 100 nL/min using a nanoflow liquid chromatography system.
##'   The effective gradient was 70 mins, allowing 16 cells per day.
##' - **Ionization**: Peptides were analyzed using an Orbitrap Eclipse mass
##'   spectrometer with a FAIMS Pro interface. FAIMS compensation voltages of
##'   -55 and -70 V were applied, with a 1-second cycle time for both voltages.
##' - **Mass spectrometry**: MS spectra were acquired with the Orbitrap
##'   analyzer, while MS/MS spectra were acquired with a linear ion trap
##'   analyzer. The maximum ion injection time for MS/MS was 200 ms.
##' - **Data analysis**: MS raw files were searched against the UniProt
##'   human protein database and an in-house contamination database
##'   using Proteome Discoverer(2.4). Label-free quantification was based on
##'   peak intensity with the match-between-runs (MBR) feature enabled.
##'
##' @section Data collection:
##'
##' The oocyte protein data shared by the author and it is accessible from the
##' [Shared File](https://biopic-my.sharepoint.cn/:x:/g/personal/humo_biopic_pku_edu_cn/EfX4CHedVopLuSx2OJNj6LABdESGNdKz4Eh8Zawvd-fNNQ?e=E5m09k&xsdata=MDV8MDJ8ZW5lcy5heWFyQHVjbG91dmFpbi5iZXxjYjY2M2MwYzNjMDY0YjZhNjc1NTA4ZGM4YzMzNjc1YXw3YWIwOTBkNGZhMmU0ZWNmYmM3YzQxMjdiNGQ1ODJlY3wxfDB8NjM4NTM5Mzk5NjI1Mzg1NDQ3fFVua25vd258VFdGcGJHWnNiM2Q4ZXlKV0lqb2lNQzR3TGpBd01EQWlMQ0pRSWpvaVYybHVNeklpTENKQlRpSTZJazFoYVd3aUxDSlhWQ0k2TW4wPXwwfHx8&sdata=Zmt4YnZFZFViTitJRkdTc0FTK2thMjdTT0EzV2JJeS83WlZmV3R6SzdvRT0%3d)
##' The K563 protein data is accessible from the
##' GitHub (https://github.com/dionezhang/CPM/blob/master/ProteinAbundance.Rdata).
##'
##' - `DataMatrix-oocyte-20240614.csv`: normalized imputed protein matrix
##' - `ProteinAbundance.Rdata`: protein matrices (normalized, log transformed)
##'
##' We initialized an empty QFeatures object and added the corresponding
##' protein assays as [SingleCellExperiment] objects.
##'
##' The oocyte protein data were exported from the shared link as
##' (`DataMatrix-oocyte-20240614.csv`). The data were formatted to a
##' [SingleCellExperiment] object and the SampleType information were added
##' as only metadata, and stored in the `colData`. The object is then added
##' to the [QFeatures] object.
##'
##' The 562 cells protein data were downloaded from the GitHub link and loaded
##' to the memory. The `Norm` object were formatted to a [SingleCellExperiment]
##' object and the SampleType information were added as only metadata, and
##' stored in the `colData`. The object is then added to the [QFeatures] object.
##'
##' @source
##'
##' The oocyte data were downloaded from the [Shared
##'     File](https://biopic-my.sharepoint.cn/:x:/g/personal/humo_biopic_pku_edu_cn/EfX4CHedVopLuSx2OJNj6LABdESGNdKz4Eh8Zawvd-fNNQ?e=E5m09k&xsdata=MDV8MDJ8ZW5lcy5heWFyQHVjbG91dmFpbi5iZXxjYjY2M2MwYzNjMDY0YjZhNjc1NTA4ZGM4YzMzNjc1YXw3YWIwOTBkNGZhMmU0ZWNmYmM3YzQxMjdiNGQ1ODJlY3wxfDB8NjM4NTM5Mzk5NjI1Mzg1NDQ3fFVua25vd258VFdGcGJHWnNiM2Q4ZXlKV0lqb2lNQzR3TGpBd01EQWlMQ0pRSWpvaVYybHVNeklpTENKQlRpSTZJazFoYVd3aUxDSlhWQ0k2TW4wPXwwfHx8&sdata=Zmt4YnZFZFViTitJRkdTc0FTK2thMjdTT0EzV2JJeS83WlZmV3R6SzdvRT0%3d)
##'     The K563 cells protein data downloaded from the GitHub
##'     (https://github.com/dionezhang/CPM/blob/master/ProteinAbundance.Rdata)
##'     The raw data and the quantification data can also be found in the
##'     MassIVE repository `MSV000089625`: ftp://MSV000089625@massive.ucsd.edu/.
##'
##' @references
##'
##' - Hu, M., Zhang, Y., Yuan, Y., Ma, W., Zheng, Y., Gu, Q., & Xie, X. S. 2023.
##'   “Correlated protein modules revealing functional coordination of
##'   interacting proteins are detected by single-cell proteomics.”. The Journal
##'   of Physical Chemistry B, ([link to
##'   article](https://doi.org/10.1021/acs.jpcb.3c00014)).
##'
##' @aliases hu2023
##' @aliases hu2023_K562
##' @aliases hu2023_oocyte
##'
##' @examples
##' \donttest{
##' hu2023_oocyte()
##' hu2023_K562()
##' }
##'
##' @keywords datasets
##'
"hu2023_K562"