# Woo et al. 2022 (Cell Syst.): 26 primary human lung cells

Single-cell proteomics data from dissociated primary human lung cells.
The data were acquired using the TIFF (transfer identification based on
FAIMS filtering) acquisition method. The data contain 26 single cells.

## Usage

``` r
woo2022_lung
```

## Format

A
[QFeatures::QFeatures](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.html)
object with 5 assays, each assay being a
[SingleCellExperiment::SingleCellExperiment](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
object:

- `peptides_[intensity or LFQ]`: 2 assays containing peptide quantities
  or normalized quantities using the maxLFQ method as computed by
  MaxQuant.

- `proteins_[intensity or iBAQ or LFQ]`: 3 assays containing protein
  quantities or normalized proteins using the iBAQ or maxLFQ methods as
  computed by MaxQuant.

Sample annotation is stored in `colData(woo_lung())`.

## Source

The peptide and protein data can be downloaded from the MASSIVE
repository MSV000085937

## Acquisition protocol

The data were acquired using the following setup. More information can
be found in the source article (see `References`).

- **Sample isolation**: primary human lung cells were dissociated
  following the protocol in Bandyopadhyay et al., 2018. The cells were
  sorted using the Influx II cell sorter and deposited on a nanoPOTS
  chip.

- **Sample preparation**: cells are lysed using using a DDM+DTT lysis
  and reduction buffer. The proteins are alkylated with IAA and digested
  with LysC and trypsin. Samples are then acidified with FA, vacuum
  dried and stored in freezer until data acquisition.

- **Liquid chromatography**: peptides are loaded using an in-house
  autosampler (Williams et al. 2020). The samples are concentrated
  through a SPE column (4cm x 100µm i.d. packed with 5µm C18) with
  microflow LC pump. The peptides are then eluted from an LC column
  (25cm x 50 µm i.d. packed with 1.7µm C18) from a 60 min gradient
  (100nL/min).

- **Mass spectrometry**: MS/MS was performed on an Orbitrap Fusion Lumos
  Tribrid MS with FAIMSpro coupled to a 2.4 kV ESI. FAIMS setup: 4-CV
  method (-45, -55, -65, -75 V). MS1 setup: resolution = 120.000, range
  = 350-1500 m/z,AGC target of 1E6, accumulation of 254ms. MS2 setup:
  30% HCD, resolution AGC 2E4, accumulation of 254ms.

- **Raw data processing**: preprocessing using Maxquant v1.6.2.10 that
  use Andromeda search engine (with UniProtKB 2016-21-29). MBR was
  enabled.

## Data collection

All data were collected from the MASSIVE repository (accession ID:
MSV000085937).

The peptide and protein data were extracted from the
`peptides_nondepleted_Lung_scProteomics.txt` or
`proteinGroups_nondepleted_Lung_scProteomics.txt` files, respectively,
in the `NonDepleted_Lung_SingleCellProteomics` folders.

The tables were split so that intensities, maxLFQ, and iBAQ data are
contained in separate tables. Tables are then converted to
[SingleCellExperiment::SingleCellExperiment](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
objects. Sample annotations were inferred from the sample names. All
data is combined in a
[QFeatures::QFeatures](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.html)
object.
[QFeatures::AssayLinks](https://rformassspectrometry.github.io/QFeatures/reference/AssayLinks.html)
were stored between peptide assays and their corresponding proteins
assays based on the leading razor protein (hence only unique peptides
are linked to proteins).

The script to reproduce the `QFeatures` object is available at
`system.file("scripts", "make-data_woo2022_lung.R", package = "scpdata")`

## References

**Source article**: Woo, Jongmin, Geremy C. Clair, Sarah M. Williams,
Song Feng, Chia-Feng Tsai, Ronald J. Moore, William B. Chrisler, et al.
2022. “Three-Dimensional Feature Matching Improves Coverage for
Single-Cell Proteomics Based on Ion Mobility Filtering.” Cell Systems 13
(5): 426–34.e4. ([link to
article](http://dx.doi.org/10.1016/j.cels.2022.02.003)).

## Examples

``` r
# \donttest{
woo2022_lung()
#> see ?scpdata and browseVignettes('scpdata') for documentation
#> loading from cache
#> An instance of class QFeatures (type: scp) with 5 sets:
#> 
#>  [1] peptides_intensity: SingleCellExperiment with 5415 rows and 26 columns 
#>  [2] peptides_LFQ: SingleCellExperiment with 5415 rows and 26 columns 
#>  [3] proteins_intensity: SingleCellExperiment with 1243 rows and 26 columns 
#>  [4] proteins_iBAQ: SingleCellExperiment with 1243 rows and 26 columns 
#>  [5] proteins_LFQ: SingleCellExperiment with 1243 rows and 26 columns 
# }
```
