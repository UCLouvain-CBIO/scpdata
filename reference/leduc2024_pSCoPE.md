# Leduc et al. 2024 - pSCoPE (Nat. Protoc.): monocyte, melanoma and PDAC cells

Single cell proteomics data acquired by the Slavov Lab. This is an
update of `leduc2022_pSCoPE`. The dataset contains quantitative
information from three cell lines (U-937 monocyte, WM-989 Melanoma and
CPAFII PDAC) at precursor, peptide and protein levels. This version of
the data was acquired using the pSCoPE MS acquisition approach.

## Format

A
[QFeatures](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.html)
object with 64 sets, each set being a
[SingleCellExperiment](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
object:

- Sets 1-61: log-transformed precursor-level single-cell data acquired
  with a 32-plex TMTpro protocol, hence these assays contain 32 columns.
  Columns contain quantitative information from single-cell channels,
  carrier channels, reference channels, empty (negative control)
  channels and unused channels.

- `peptides`: peptide-level data containing quantitative data for 13427
  peptides in 1769 single cells.

- `proteins`: protein-level data aggregated based on
  `Leading.razor.protein`, containing quantitative data for 1909
  proteins in 1769 single cells.

- `imported_proteins`: protein-level table imported from the authors'
  protein report, containing 1775 proteins in 1526 single cells.

The `colData(leduc2024_pSCoPE())` contains cell type annotation and
batch annotation that are common to all sets.

## Source

The links to the data were found on the [Slavov Lab
website](https://scp.slavovlab.net/Leduc_et_al_2023). The data were
downloaded from the [Google drive
folder](https://drive.google.com/drive/folders/10LYxtc7qNdO06z8FKtf0L9UB4IAta68h).
The raw data and quantification data can also be found in the MassIVE
repository `MSV000094207`: ftp://massive.ucsd.edu/MSV000094207.

## Acquisition protocol

The data were acquired using the following setup. More information can
be found in the source article (see `References`).

- **Cell isolation**: CellenONE cell sorting.

- **Sample preparation** performed using the improved SCoPE2 protocol
  using the CellenONE liquid handling system: nPOP cell lysis (DMSO) +
  trypsin digestion + TMT-32plex labelling and pooling. The DDA data
  were acquired following the procedure outlined in Huffman et al.
  (2022). **Separation**: online nLC (Neo Vanquish HPLC) with a 25cm x
  75um IonOpticks Aurora column (active gradient ~24 min).

- **Mass spectrometry**: Orbitrap Exploris 480 mass spectrometer.

- **Data analysis**: MaxQuant.

## Data collection

The PSM data were collected from a massIVE repository (see `Source`
section). The archive `f.MSV000094207/search/All_1p_FDR.zip` contains
the following file of interest:

- `evidence.txt`: PSM-level quantitative data

Sample metadata were collected from a shared Google Drive folder that is
accessible from the SlavovLab website (see `Source` section). The folder
contains the following files of interest:

- `TMT_SingleCell_Meta.csv`: sample-level annotations

- `TMT_ProteinxSingleCell.csv`: processed protein-level data table
  provided by the authors, containing the `imported_proteins` data

We formatted sample-level annotations to ensure columns match with those
of the quantification tables. Sample annotations containing `NA` values
for `neg` samples were manually filled in based on the `ID` column in
order to circumvent duplicate sample naming errors. These changes are
incorporated in the file `TMT_ProteinxSingleCell_modified.csv`,
available on [Zenodo](https://doi.org/10.5281/zenodo.24020112).

Both quantification and annotation tables were then combined in a single
[QFeatures](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.html)
object using
[`QFeatures::readQFeatures()`](https://rformassspectrometry.github.io/QFeatures/reference/readQFeatures.html).

The script to reproduce the QFeatures object is available at
inst/scripts/make-data_leduc2024_pSCoPE.R.

## References

Andrew Leduc, Luke Khoury, Joshua Cantlon, Saad Khan & Nikolai Slavov.
2024. “Massively parallel sample preparation for multiplexed single-cell
proteomics using nPOP.” Nature Protocols 19: 3750–3776. [Link to
article](https://doi.org/10.1038/s41596-024-01033-8).

Andrew Leduc, Gray Huffman, Joshua Cantlon, Saad Khan & Nikolai Slavov,
2022. “Exploring functional protein covariation across single cells
using nPOP.” Genome Biology 23, 261. [Link to
article](https://doi.org/10.1186/s13059-022-02817-5)

Gray Huffman, Andrew Leduc, Christoph Wichmann, Marco di Gioia,
Francesco Borriello, Harrison Specht, Jason Derks, et al. 2023.
"Prioritized mass spectrometry increases the depth, sensitivity and data
completeness of single-cell proteomics." Nature Methods 20: 714-722.
[Link to article](https://doi.org/10.1038/s41592-023-01830-1)

## See also

[leduc2024_plexDIA](https://UCLouvain-CBIO.github.io/scpdata/reference/leduc2024_plexDIA.md)

## Examples

``` r
# \donttest{
leduc2024_pSCoPE()
#> Error in leduc2024_pSCoPE(): "leduc2024_pSCoPE" not found in ExperimentHub
# }
```
