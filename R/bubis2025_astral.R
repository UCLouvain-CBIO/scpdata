##' Bubis et al. 2023 (Nat Methods): Challenging the Astral mass analyzer.
##'
##' Single-cell (or single-cell range peptide dilution) proteomics data
##' acquired on the Thermo Scientific Orbitrap Astral mass spectrometer.
##' The dataset contains single-cell range peptide dilution from Hela or
##' K562 cell lines, single-cell data from A549 and H460 cell lines
##' (epithelial-like human lung cancer cells), as well as trophectoderm-
##' like (TE-like) and naive human pluripotent stem (hPS) cells. The
##' dataset contains quantitative information at the PSM, peptide and
##' protein levels. The raw data were analysed using Spectronaut version
##' 18.
##'
##' @format A [QFeatures] object with 602 sets, each set being a
##' [SingleCellExperiment] object:
##'
##' - Sets 1-594: PSM data from the Spectronaut peptides report files
##'   with LFQ quantities from the `PEP.MS1Quantity` column.
##' - Sets 595-598: Peptide data resulting from the PSM to peptide
##'   aggregation of the 594 PSM set. Resulting peptide sets were
##'   joined into assay sets based on the figure they were used for (1,
##'   2, 4 or 5).
##' - Sets 599-602: Protein data from the Spectronaut protein report
##'   files with LFQ quantities from the `PG.Quantity` column.
##'
##' The `colData(bubis2025())` contains cell type annotation, search
##' strategy, LC method throughput, FAIMS compensation voltage and figure
##' in which the sample was used. The description of the `rowData` fields
##' can be found in the [`Spectronaut` user manual](https://biognosys.com/content/uploads/2025/06/Spectronaut-20-Manual.pdf).
##'
##' @section Acquisition protocol:
##'
##' The data were acquired using the following setup. More information
##' can be found in the source article (see *References*).
##'
##' - **Cell isolation**: A549 and H460 cells were isolated into 384-well
##'   plates using the cellenONE robot. TE-like cells and hSP cells were
##'   sorted into 384-well plates using the FACS Aria III (BD) controlled
##'   by the DiVa software version 9.0.1. TE-like cells were selected
##'   based on their TROP2 positivity while naive hPS cells were selected
##'   based on their SUSD2 positivity. Whether cells were sorted using
##'   cellenONE or FACS, the wells were pre-filled with 1 microL of
##'   lysis/digestion buffer.
##' - **Sample preparation** Lysis and digestion were performed in the
##'   plates for 30 min at 50°C, then stopped by adding 0.1% TFA. Liquid
##'   dispensing was performed using the cellenONE robot. The plates were
##'   stored at -20°C until injection.
##' - **Liquid chromatography**: Samples were injected directly from the
##'   384-well plate. Chromatographic separation of peptides was
##'   performed on a vanquish Neo UHPLC system connected to an Aurora
##'   Ultimate TS 25-cm nanoflow UHPLC column with an integrated emitter
##'   (Ion Optics) at 50 °C in direct injection mode. Peptide separation
##'   was performed at 30-80 SPD, with details provided in the `colData`
##'   for each sample.
##' - **Mass spectrometry**: Acquisition was conducted using an Orbitrap
##'   Astral mass spectrometer equipped with the FAIMS Pro or FAIMS Pro
##'   Duo interface and an EASY-Spray source. FAIMS compensation voltage
##'   (CV) was set to -48 V except for the samples used in fig. 2, for
##'   which different CV were tested. The exact CV used for each sample
##'   is specified in the `colData`. MS1 spectra were recorded using the
##'   Orbitrap analyser at a resolution of 240,000, from m/z 400 to 900,
##'   using an automated gain control (AGC) target of 500 % and a maximum
##'   injection time of 100 ms. For MS2 in DIA mode using the Astral
##'   analyzer, non-overlapping isolation windows of 20 m/z. A scan range
##'   of m/z 400 to 800 was selected. Precursor accumulation time was set
##'   to 60 ms and the AGC target to 800%.
##' - **Raw data processing**: Raw files were processed using
##'   Spectronaut version 18.6.231227.55695. For results marked as
##'   ‘method evaluation’, DirectDIA+ was used in method evaluation mode
##'   without cross-normalisation. ‘DirectDIA+’ indicates that results
##'   were analyzed in DirectDIA+ mode. For library searches, higher-
##'   input DIA recorded results were first analyzed in DirectDIA+ mode,
##'   and a library was created from those files using Spectronaut.
##'   Quantification was performed at the MS1 level. Carbamidomethylation
##'   of cysteines as a static modification was removed for single-cell
##'   searches. Factory settings were used for all analyses and for
##'   library generation unless otherwise indicated. ‘Optimized settings’
##'   (OS) indicates that searches were performed with a more stringent
##'   cutoff value of 0.01. The search strategy used for each sample is
##'   specified in the `colData`. Searches were performed against the
##'   human proteome (UniProt proteome UP000005640) and the CRAPome.
##'
##' @section Data collection:
##'
##' The MS raw proteomics data, Spectronaut search results and fasta
##' files used are available on the ProteomeXchange Consortium via the
##' PRIDE partner repository with the dataset identifier [PXD049412](https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD049412).
##' Only files related to single-cells or single-cell range peptide
##' dilutions were used for this dataset (fig. 1, 2, 4 and 5). Sample
##' annotations were manually generated using the file names and
##' descriptions from the source article.
##'
##' PSM and proteins data were extracted from the
##' `_Report_Peptides_JB_Pivot.tsv` and `_Report_Protein_JB_Pivot.tsv`
##' files respectively.
##'
##' Sample annotations and PSM features tables were combined in a single
##' [QFeatures] object using the [scp::readSCP()] function. Peptide data
##' were obtained by aggregating PSM data to the peptide level. Resulting
##' peptides sets were joined according to the figure for which they were
##' used in the source article. Individual peptides sets were discarded.
##' Protein data were made in a separate [QFeatures] object, and sets
##' were joined based on the figure they are related to. Combined sets
##' were added to the first [QFeatures] object containing PSM and peptide
##' data, and links between peptide and protein sets were set manually.
##'
##' Note that the [QFeatures] object has not been further processed and
##' has therefore not been normalized, log-transformed or batch-
##' corrected.
##'
##' @source The PSM and protein data, as well as sample annotations, can
##'     be downloaded on the ProteomeXchange Consortium via the PRIDE
##'     partner repository with the dataset identifier [PXD049412](https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD049412).
##'
##' @references
##'
##' Julia A. Bubis, Tabiwang N. Arrey, Eugen Damoc, Bernard Delanghe,
##' Jana Slovakova, Theresa M. Sommer, Harunobu Kagawa, Peter Pichler,
##' Nicolas Rivron, Karl Mechtler and Manuel Matzinger.
##'
##' 2025. "Challenging the Astral mass analyzer to quantify up to 5,300 proteins per single cell at unseen accuracy to uncover cellular heterogeneity"
##' Nature Methods  https://doi.org/10.1038/s41592-024-02559-1
##' DOI:[10.1038/s41592-024-02559-1](https://doi.org/10.1038/s41592-024-02559-1)
##'
##' @examples
##' \donttest{
##' bubis2025()
##' }
##'
##' @keywords datasets
##'
"bubis2025"
