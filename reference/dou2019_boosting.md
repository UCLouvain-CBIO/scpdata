# Dou et al. 2019 (Anal. Chem.): testing boosting ratios

Single-cell proteomics using nanoPOTS combined with TMT isobaric
labeling. It contains quantitative information at PSM and protein level.
The cell types are either "Raw" (macrophage cells), "C10" (epihelial
cells), or "SVEC" (endothelial cells). Each cell is replicated 2 or 3
times. Each cell type was run using 3 levels of boosting: 0 ng (no
boosting), 5 ng or 50 ng. When boosting was applied, 1 reference well
and 1 boosting well were added, otherwise 1 empty well was added. Each
boosting setting (0ng, 5ng, 50ng) was run in duplicate.

## Usage

``` r
dou2019_boosting
```

## Format

A
[QFeatures::QFeatures](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.html)
object with 7 assays, each assay being a
[SingleCellExperiment::SingleCellExperiment](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
object:

- `Boosting_X_run_Y`: PSM data with 10 columns corresponding to the
  TMT-10plex channels. The `X` indicates the boosting amount (0ng, 5ng
  or 50ng) and `Y` indicates the run number (1 or 2).

- `peptides`: peptide data containing quantitative data for 13,462
  peptides in 60 samples (run 1 and run 2 combined).

- `proteins`: protein data containing quantitative data for 1436
  proteins and 60 samples (all runs combined).

Sample annotation is stored in `colData(dou2019_boosting())`.

## Source

The PSM data can be downloaded from the massIVE repository MSV000084110.
FTP link: ftp://massive.ucsd.edu/MSV000084110/

The protein data can be downloaded from the [ACS
Publications](https://pubs.acs.org/doi/10.1021/acs.analchem.9b03349)
website (Supporting information section).

## Acquisition protocol

The data were acquired using the following setup. More information can
be found in the source article (see `References`).

- **Cell isolation**: single-cells from the three murine cell lines were
  isolated using FACS (BD Influx II cell sorter ). Boosting sample were
  prepared (presumably in bulk) from 1:1:1 mix of the three cell lines.

- **Sample preparation** performed using the nanoPOTs device. Protein
  extraction (DMM + TCEAP) + alkylation (IAA) + Lys-C digestion +
  trypsin digestion + TMT-10plex labeling and pooling.

- **Separation**: nanoLC (Dionex UltiMate with an in-house packed 50cm x
  30um LC columns; 50nL/min)

- **Ionization**: ESI (2,000V)

- **Mass spectrometry**: Thermo Fisher Orbitrap Fusion Lumos Tribrid
  (MS1 accumulation time = 50ms; MS1 resolution = 120,000; MS1 AGC =
  1E6; MS2 accumulation time = 246ms; MS2 resolution = 60,000; MS2 AGC =
  1E5)

- **Data analysis**: MS-GF+ + MASIC (v3.0.7111) + RomicsProcessor
  (custom R package)

## Data collection

The PSM data were collected from the MassIVE repository MSV000084110
(see `Source` section). The downloaded files are:

- `Boosting_*ng_run_*_msgfplus.mzid`: the MS-GF+ identification result
  files.

- `Boosting_*ng_run_*_ReporterIons.txt`: the MASIC quantification result
  files.

For each batch, the quantification and identification data were combined
based on the scan number (common to both data sets). The combined
datasets for the different runs were then concatenated feature-wise. To
avoid data duplication due to ambiguous matching of spectra to peptides
or ambiguous mapping of peptides to proteins, we combined ambiguous
peptides to peptides groups and proteins to protein groups. Feature
annotations that are not common within a peptide or protein group are
are separated by a `;`. The sample annotation table was manually created
based on the available information provided in the article. The data
were then converted to a
[QFeatures::QFeatures](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.html)
object using the
[`scp::readSCP()`](https://UCLouvain-CBIO.github.io/scp/reference/readSCP.html)
function.

We generated the peptide data. First, we removed PSM matched to
contaminants or decoy peptides and ensured a 1% FDR. We aggregated the
PSM to peptides based on the peptide (or peptide group) sequence(s)
using the median PSM instenity. The peptide data for the different runs
were then joined in a single assay (see
[QFeatures::joinAssays](https://rformassspectrometry.github.io/QFeatures/reference/joinAssays.html)),
again based on the peptide sequence(s). We then removed the peptide
groups. Links between the peptide and the PSM data were created using
[QFeatures::addAssayLink](https://rformassspectrometry.github.io/QFeatures/reference/AssayLinks.html).
Note that links between PSM and peptide groups are not stored.

The protein data were downloaded from `Supporting information` section
from the publisher's website (see `Sources`). The data is supplied as an
Excel file `ac9b03349_si_004.xlsx`. The file contains 7 sheets from
which we took the 2nd, 4th and 6th sheets (named
`01 - No Boost raw data`, `03 - 5ng boost raw data`,
`05 - 50ng boost raw data`, respectively). The sheets contain the
combined protein data for the duplicate runs given the boosting amount.
We joined the data for all boosting ration based on the protein name and
converted the data to a
[SingleCellExperiment::SingleCellExperiment](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
object. We then added the object as a new assay in the
[QFeatures::QFeatures](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.html)
dataset (containing the PSM data). Links between the proteins and the
corresponding PSM were created. Note that links to protein groups are
not stored.

## References

Dou, Maowei, Geremy Clair, Chia-Feng Tsai, Kerui Xu, William B.
Chrisler, Ryan L. Sontag, Rui Zhao, et al. 2019. “High-Throughput Single
Cell Proteomics Enabled by Multiplex Isobaric Labeling in a Nanodroplet
Sample Preparation Platform.” Analytical Chemistry, September ([link to
article](https://doi.org/10.1021/acs.analchem.9b03349)).

## See also

[dou2019_lysates](https://UCLouvain-CBIO.github.io/scpdata/reference/dou2019_lysates.md),
[dou2019_mouse](https://UCLouvain-CBIO.github.io/scpdata/reference/dou2019_mouse.md)

## Examples

``` r
# \donttest{
dou2019_boosting()
#> see ?scpdata and browseVignettes('scpdata') for documentation
#> loading from cache
#> An instance of class QFeatures (type: scp) with 8 sets:
#> 
#>  [1] Boosting_0ng_run_1: SingleCellExperiment with 12258 rows and 10 columns 
#>  [2] Boosting_0ng_run_2: SingleCellExperiment with 12611 rows and 10 columns 
#>  [3] Boosting_50ng_run_1: SingleCellExperiment with 27746 rows and 10 columns 
#>  ...
#>  [6] Boosting_5ng_run_2: SingleCellExperiment with 20435 rows and 10 columns 
#>  [7] peptides: SingleCellExperiment with 13462 rows and 60 columns 
#>  [8] proteins: SingleCellExperiment with 1436 rows and 60 columns 
# }
```
