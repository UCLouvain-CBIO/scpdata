# Leduc et al. 2024 - plexDIA (Nat. Protoc.): monocyte, melanoma and PDAC cells

Single cell proteomics data acquired by the Slavov Lab. This is an
update of `leduc2022_plexDIA`. The dataset contains quantitative
information from three cell lines (U-937 monocyte, WM-989 Melanoma and
CPAFII PDAC) at precursor, peptide and protein levels. This version of
the data was acquired using the plexDIA MS acquisition protocol.

## Format

A
[QFeatures](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.html)
object with 50 sets, each set being a
[SingleCellExperiment](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
object:

- Sets 1-47: log-transformed precursor-level single-cell data acquired
  using a 3-plex mTRAQ protocol, hence these assays contain 3 columns.
  Columns contain quantitative information from single cells or negative
  control samples.

- `peptides`: peptide-level data containing quantitative data for 94080
  peptides in 141 single cells.

- `proteins`: protein-level data aggregated based on `Protein.Ids`,
  containing quantitative data for 8046 proteins in 141 single cells.

- `imported_proteins`: protein-level table imported from the authors'
  protein report, containing 4619 proteins in 121 single cells.

The `colData(leduc2024_plex_DIA())` contains cell type annotation and
batch annotation that are common to all sets.

## Source

The links to the data were found on the [Slavov Lab
website](https://scp.slavovlab.net/Leduc_et_al_2023). The data were
downloaded from the [Google drive
folder](https://drive.google.com/drive/folders/1gxBS7PlWSGzmql301hYglVgTVq1YRz-z).
The raw data and quantification data can also be found in the MassIVE
repository `MSV000093494`: ftp://massive.ucsd.edu/MSV000093494.

## Acquisition protocol

The data were acquired using the following setup. More information can
be found in the source article (see `References`).

- **Cell isolation**: CellenONE cell sorting.

- **Sample preparation**: performed using the improved SCoPE2 protocol
  using the CellenONE liquid handling system: nPOP cell lysis (DMSO) +
  trypsin digestion + mTRAQ 3-plex labelling.

- **Separation**: online nLC (NanoElute2) with a 25cm x 75um IonOpticks
  Aurora column and captive spray fitting (LC gradient 4-38% buffer B).

- **Mass spectrometry**: timsTOF Ultra mass spectrometer.

- **Data analysis**: DIA-NN.

## Data collection

The PSM data were collected from a massIVE repository (see `Source`
section). The folder `f.MSV000093494/search` contains the following file
of interest:

- `report.tsv`: PSM-level quantitative data

Sample metadata were collected from a shared Google Drive folder that is
accessible from the SlavovLab website (see `Source` section). The folder
contains the following files of interest:

- `plexDIA_SingleCell_Meta.csv`: sample-level annotations

- `plexDIA_SingleCellxProtein.csv`: processed protein-level data table
  provided by the authors, containing the `imported_proteins` data

We formatted sample-level annotations to ensure columns match with those
of the quantification tables. Sample annotations containing `NA` values
for `neg` samples were manually filled in based on the `ID` column in
order to circumvent duplicate sample naming errors. These changes are
incorporated in the file `plexDIA_SingleCellxProtein_modified.csv`,
available on [Zenodo](https://doi.org/0.5281/zenodo.24020222).

Both quantification and annotation tables were then combined in a single
[QFeatures](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.html)
object using
[`QFeatures::readQFeatures()`](https://rformassspectrometry.github.io/QFeatures/reference/readQFeatures.html).

The script to reproduce the QFeatures object is available at
inst/scripts/make-data_leduc2024_plexDIA.R.

## References

Andrew Leduc, Luke Khoury, Joshua Cantlon, Saad Khan & Nikolai Slavov.
2024. “Massively parallel sample preparation for multiplexed single-cell
proteomics using nPOP.” Nature Protocols 19: 3750–3776. [Link to
article](https://doi.org/10.1038/s41596-024-01033-8).

Andrew Leduc, Gray Huffman, Joshua Cantlon, Saad Khan & Nikolai Slavov,
2022. “Exploring functional protein covariation across single cells
using nPOP.” Genome Biology 23, 261. [Link to
article](https://doi.org/10.1186/s13059-022-02817-5)

## See also

[leduc2024_pSCoPE](https://UCLouvain-CBIO.github.io/scpdata/reference/leduc2024_pSCoPE.md)

## Examples

``` r
# \donttest{
leduc2024_plexDIA()
#> Error in leduc2024_plexDIA(): "leduc2024_plexDIA" not found in ExperimentHub
# }
```
