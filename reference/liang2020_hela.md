# Liang et al. 2020 (Anal. Chem.): HeLa cells (MaxQuant preprocessing)

Single-cell proteomics data from HeLa cells using the autoPOTS
acquisition workflow. The samples contain either no cells (blanks), 1
cell, 10 cells, 150 cells or 500 cells. Samples containing between 0 and
10 cells are isolated using micro-pipetting while samples containing
between 150 and 500 cells were prepared using dilution of a bulk sample.

## Usage

``` r
liang2020_hela
```

## Format

A
[QFeatures::QFeatures](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.html)
object with 17 assays, each assay being a
[SingleCellExperiment::SingleCellExperiment](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
object:

- `HeLa_*`: 15 assays containing PSM data.

- `peptides`: quantitative data for 48705 peptides in 15 samples (all
  runs are combined).

- `proteins`: quantitative data for 3970 protein groups in 15 samples
  (all runs combined).

Sample annotation is stored in `colData(liang2020_hela())`.

## Source

The PSM data can be downloaded from the PRIDE repository PXD021882 The
source link is:
http://ftp.pride.ebi.ac.uk/pride/data/archive/2020/12/PXD021882/

## Acquisition protocol

The data were acquired using the following setup. More information can
be found in the source article (see `References`).

- **Cell isolation**: The HeLa cells come from a commercially available
  cell line. Samples containing between 0 and 10 cells were isolated
  using micro-manipulation and the counts were validated using a
  microscope. Samples containing between 150 and 500 cells were prepared
  by diluting a bulk sample and the exact counts were evaluated by
  obtaining phtotmicrographs.

- **Sample preparation** performed using the autoPOTS worflow that
  relied on the OT-2 pipeting robot. Cell are lysed using sonication.
  Samples are then processed by successive incubation with DTT
  (reduction), then IAA (alkylation), then Lys-C and trypsin (protein
  digestion).

- **Separation**: Samples were injected on the column using a modified
  Ultimate WPS-3000 TPL autosampler coupled to an UltiMate 3000 RSLCnano
  pump. The LC column is a home-packed nanoLC column (45cm x 30um;
  40nL/min)

- **Ionization**: Nanospray Flex ion source (2,000V)

- **Mass spectrometry**: Orbitrap Exploris 480. MS1 settings:
  accumulation time = 250 ms (0-10 cells) or 100 ms (150-500 cells);
  resolution = 120,000; AGC = 100\\ duration = 90 s (0-10 cells) or 60 s
  (150-500 cells) ; accumulation time = 500 ms (0-1 cell), 250 ms (10
  cells), 100 ms (150 cells) or 50 ms (500 cells); resolution = 60,000
  (0-10 cells) or 30,000 (150-500 cells); AGC = 5E3 (0-1 cells) or 1E4
  (10-500 cells).

- **Data analysis**: MaxQuant (v1.6.7.0) and the search database is
  Swiss-Prot (July 2020).

## Data collection

All data were collected from the PRIDE repository (accession ID:
PXD021882).

The sample annotations were collected from the methods section and from
table S3 in the paper.

The PSM data were found in the `evidence.txt` file. The data were
converted to a
[QFeatures::QFeatures](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.html)
object using the
[`scp::readSCP()`](https://UCLouvain-CBIO.github.io/scp/reference/readSCP.html)
function.

The peptide data were found in the `peptides.txt` file. The column names
holding the quantitative data were adapted to match the sample names in
the
[QFeatures::QFeatures](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.html)
object. The data were then converted to a
[SingleCellExperiment::SingleCellExperiment](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
object and then inserted in the
[QFeatures::QFeatures](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.html)
object. Links between the PSMs and the peptides were added

A similar procedure was applied to the protein data. The data were found
in the `proteinGroups.txt` file. The column names were adapted, the data
were converted to a
[SingleCellExperiment::SingleCellExperiment](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
object and then inserted in the
[QFeatures::QFeatures](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.html)
object. Links between the peptides and the proteins were added

## References

Liang, Yiran, Hayden Acor, Michaela A. McCown, Andikan J. Nwosu, Hannah
Boekweg, Nathaniel B. Axtell, Thy Truong, Yongzheng Cong, Samuel H.
Payne, and Ryan T. Kelly. 2020. “Fully Automated Sample Processing and
Analysis Workflow for Low-Input Proteome Profiling.” Analytical
Chemistry, December. ([link to
article](https://doi.org/10.1021/acs.analchem.0c04240)).

## Examples

``` r
# \donttest{
liang2020_hela()
#> see ?scpdata and browseVignettes('scpdata') for documentation
#> loading from cache
#> An instance of class QFeatures (type: scp) with 17 sets:
#> 
#>  [1] HeLa_0cell_B16: SingleCellExperiment with 1153 rows and 1 columns 
#>  [2] HeLa_0cell_B18: SingleCellExperiment with 1214 rows and 1 columns 
#>  [3] HeLa_0cell_B4: SingleCellExperiment with 1366 rows and 1 columns 
#>  ...
#>  [15] HeLa_500cell_K8: SingleCellExperiment with 35054 rows and 1 columns 
#>  [16] peptides: SingleCellExperiment with 48705 rows and 15 columns 
#>  [17] proteins: SingleCellExperiment with 3970 rows and 15 columns 
# }
```
