# Zhu et al. 2018 (Mol. Cel. Prot.): rat brain laser dissections

Near single-cell proteomics data of laser captured micro-dissection
samples. The samples are 24 brain sections from rat pups (day 17). The
slices are 12 um thick squares of either 50, 100, or 200 um width. 5
samples were dissected from the corpus callum (`CC`), 4 samples were
dissected from the corpus collosum (`CP`), 13 samples were extracted
from the cerebral cortex (`CTX`), and 2 samples are labeled as (`Mix`).

## Usage

``` r
zhu2018MCP
```

## Format

A
[QFeatures::QFeatures](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.html)
object with 4 assays, each assay being a
[SingleCellExperiment::SingleCellExperiment](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
object:

- `peptides`: quantitative information for 13,055 peptides from 24
  samples

- `proteins_intensity`: protein intensities for 2,257 proteins from 24
  samples

- `proteins_LFQ`: LFQ intensities for 2,257 proteins from 24 samples

- `proteins_iBAQ`: iBAQ values for 2,257 proteins from 24 samples

Sample annotation is stored in `colData(zhu2018MCP())`.

## Source

The PSM data can be downloaded from the PRIDE repository PXD008844. FTP
link ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2018/07/PXD008844

## Acquisition protocol

The data were acquired using the following setup. More information can
be found in the original article (see `References`).

- **Cell isolation**: brain patches were collected using laser-capture
  microdissection (PALM MicroBeam) on flash frozen rat (*Rattus
  norvergicus*) brain tissues. Note that the samples were stained with
  H&E before dissection for histological analysis. DMSO is used as
  sample collection solution

- **Sample preparation** performed using the nanoPOTs device: DMSO
  evaporation + protein extraction (DMM + DTT) + alkylation (IAA)

  - Lys-C digestion + trypsin digestion.

- **Separation**: nanoLC (Dionex UltiMate with an in-house packed 60cm x
  30um LC columns; 50nL/min)

- **Ionization**: ESI (2,000V)

- **Mass spectrometry**: Thermo Fisher Orbitrap Fusion Lumos Tribrid
  (MS1 accumulation time = 246ms; MS1 resolution = 120,000; MS1 AGC =
  3E6). The MS/MS settings depend on the sample size, excepted for the
  AGC = 1E5. 50um (time = 502ms; resolution = 240,000), 100um (time =
  246ms; resolution = 120,000), 200um (time = 118ms; resolution =
  60,000).

- **Data analysis**: MaxQuant (v1.5.3.30) + Perseus (v1.5.6.0) + Origin
  Pro 2017

## Data collection

The data were collected from the PRIDE repository (accession ID:
PXD008844). We downloaded the `MaxQuant_Peptides.txt` and the
`MaxQuant_ProteinGroups.txt` files containing the combined
identification and quantification results. The sample annotations were
inferred from the names of columns holding the quantification data and
the information in the article. The peptides data were converted to a
[SingleCellExperiment::SingleCellExperiment](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
object. We split the protein table to separate the three types of
quantification: protein intensity, label-free quantitification (LFQ) and
intensity based absolute quantification (iBAQ). Each table is converted
to a
[SingleCellExperiment::SingleCellExperiment](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
object along with the remaining protein annotations. The 4 objects are
combined in a single
[QFeatures::QFeatures](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.html)
object and feature links are created based on the peptide leading razor
protein ID and the protein ID.

## References

Zhu, Ying, Maowei Dou, Paul D. Piehowski, Yiran Liang, Fangjun Wang,
Rosalie K. Chu, William B. Chrisler, et al. 2018. “Spatially Resolved
Proteome Mapping of Laser Capture Microdissected Tissue with Automated
Sample Transfer to Nanodroplets.” Molecular & Cellular Proteomics: MCP
17 (9): 1864–74 ([link to
article](http://dx.doi.org/10.1074/mcp.TIR118.000686)).

## Examples

``` r
# \donttest{
zhu2018MCP()
#> see ?scpdata and browseVignettes('scpdata') for documentation
#> loading from cache
#> An instance of class QFeatures (type: scp) with 4 sets:
#> 
#>  [1] peptides: SingleCellExperiment with 13055 rows and 24 columns 
#>  [2] proteins_intensity: SingleCellExperiment with 840 rows and 0 columns 
#>  [3] proteins_LFQ: data.frame with 840 rows and 0 columns 
#>  [4] proteins_iBAQ: SingleCellExperiment with 840 rows and 0 columns 
# }
```
