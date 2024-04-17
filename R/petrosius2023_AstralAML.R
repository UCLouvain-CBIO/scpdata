##' Petrosius et al. 2023 (bioRxiv): AML hierarchy on Astral.
##'
##' Single cell proteomics data from FACS sorted cells from the
##' OCI-AML8227 model. The dataset contains leukemic stem cells (LSC;
##' CD34+, CD38-), progenitor cells (CD34+, CD38+), CD38+ blasts
##' (CD34-, CD38+) and CD38- blasts (CD34-, CD38-). It contains
##' quantitative information at PSM, peptide and protein levels. Data
##' was acquired using an Orbitrap Astral mass spectrometer. Direct DIA
##' analysis was performed with Spectronaut version 17.
##'
##' @format A [QFeatures] object with 217 assays, each assay being a
##' [SingleCellExperiment] object:
##'
##' - Assays 1-215: PSM data from the Spectronaut PEPQuant file with
##'   LFQ quantities from the FG.MS1Quantity column.
##' - `peptides`: Peptide data resulting from the PSM to peptide
##'   aggregation the 215 PSM assays. Resulting peptide assays were
##'   joined into a single assay.
##' - `proteins`: Protein data from the Spectronaut PGQuant file with
##'   LFQ quantities from the PG.Quantity column.
##'
##' The `colData(petrosius2023_AstralAML())` contains cell type
##' annotation, batch annotation and FACS data. The description of the
##' `rowData` fields can be found in the
##' [`Spectronaut` user manual](https://biognosys.com/content/uploads/2023/03/Spectronaut-17_UserManual.pdf).
##'
##' @section Acquisition protocol:
##'
##' The data were acquired using the following setup. More information
##' can be found in the source article (see *References*).
##'
##' - **Cell isolation**: Cell sorting was done on a FACS Aria III or
##'   Aria II instrument, controlled by the DIVA software package and
##'   operated with a 100 microm nozzle. Cells were sorted at single-cell
##'   resolution, into a 384-well Eppendorf LoBind PCR plate containing
##'   1 microL of lysis buffer.
##' - **Sample preparation** Single-cell protein lysates were digested
##'   overnight at 37°C with 2 ng of Trypsin supplied in 1 microL of
##'   digestion buffer. Digestion was stopped by the addition of 1 microL
##'   1% (v/v) trifluoroacetic acid (TFA). All liquid dispensing was
##'   done using an I-DOT One instrument.
##' - **Liquid chromatography**: Chromatographic separation of peptides
##'   was conducted on a vanquish Neo UHPLC system connected to a 50 cm
##'   uPAC Neo Low-load and an EASY-spray. Autosampler and injection
##'   valves were configured to perform direct injections from a 384
##'   well plate using a 25 uL injection loop on 11.8 min gradients.
##' - **Mass spectrometry**: Acquisition was conducted with an Orbitrap
##'   Astral mass spectrometer operated in positive mode with the
##'   FAIMSPro interface compensation voltage set to -45 V.
##'   MS1 scans were acquired with the Orbitrap at a resolution of
##'   120,000 and a scan range of 400 to 900 m/z with normalized
##'   automatic gain control (AGC) target of 300 % and maximum
##'   injection time of 246 ms. Data independent acquisition of MS2
##'   spectra was performed in the Astral using loop control set to 0.7
##'   seconds per cycle with varying isolation window widths and
##'   injection times. Fragmentation of precursor ions was performed
##'   using higher energy collisional dissociation (HCD) using a
##'   normalized collision energy (NCE) of 25 %. AGC target was set to
##'   800 %.
##' - **Raw data processing**: Raw files were processed using
##'   Spectronaut version 17. Direct DIA analysis was performed in
##'   pipeline mode. Pulsar searches were performed without fixed
##'   modifications. N-terminal acetylation and methionine oxidation
##'   were set as variable modifications. Quantification level was set
##'   to MS1 and the quantity type set to area under the curve.
##'
##' @section Data collection:
##'
##' The data were provided by the authors and is accessible at the
##' [Dataverse](https://dataverse.uclouvain.be/dataset.xhtml?persistentId=doi:10.14428/DVN/EMAVLT)
##' The dataset ('Astral AML single-cell data from Petrosius et
##' al. 2023 preprint') contains the following files of interest:
##'
##' - `20240201_130747_PEPQuant (Normal).tsv`: the PSM level data
##' - `20240201_130747_PGQuant (Normal).tsv`: the protein level data
##' - `index_map.csv`: FACS data.
##' - `msRuns_overview.csv`: Sample annotations.
##'
##' We added the FACS data to the sample annotations in a single table.
##' Both annotations and PSM features tables are then combined in a
##' single [QFeatures] object using the [scp::readSCP()] function.
##'
##' The peptide data were obtained by aggregation of the PSM data to
##' the peptide level. All of the resulting peptides assays were joined
##' into a single assays. Individual peptides assays were discarded.
##'
##' The protein data were formatted from the `20240201_130747_PGQuant (Normal).tsv`
##' to a [SingleCellExperiment] object and the sample metadata were
##' matched to the column names and stored in the  `colData`. The
##' object is then added to the [QFeatures] object and the rows of the
##' peptide data are linked to the rows of the protein data based on
##' the protein sequence information through an `AssayLink` object.
##'
##' Note that the [QFeatures] object has not been further processed and
##' has therefore not been normalized, log-transformed or
##' batch-corrected.
##'
##' @source The PSM data, protein data and sample annotations can be
##'     downloaded from the dataset 'Astral AML single-cell data from
##'     Petrosius et al. 2023 preprint' in the
##'     [Dataverse](https://dataverse.uclouvain.be/dataset.xhtml?persistentId=doi:10.14428/DVN/EMAVLT).
##'
##' @references
##'
##' Valdemaras Petrosius, Pedro Aragon-Fernandez, Tabiwang N. Arrey,
##' Nil Üresin, Benjamin Furtwängler, Hamish Stewart, Eduard Denisov,
##' Johannes Petzoldt, Amelia C. Peterson, Christian Hock, Eugen
##' Damoc, Alexander Makarov, Vlad Zabrouskov, Bo T. Porse and Erwin
##' M. Schoof.
##' 2023. "Evaluating the capabilities of the Astral mass analyzer for single-cell proteomics."
##' biorxiv.  https://doi.org/10.1101/2023.06.06.543943
##' DOI:[10.1101/2023.06.06.543943](https://doi.org/10.1101/2023.06.06.543943)
##'
##' @examples
##' \donttest{
##' petrosius2023_AstralAML()
##' }
##'
##' @keywords datasets
##'
"petrosius2023_AstralAML"
