# Single-Cell Proteomics Data Package

The `scpdata` package distributes mass spectrometry-based single-cell
proteomics datasets. The datasets were collected from published work and
formatted to a standardized data framework. The `scp` frameworks stores
the expression data for different MS levels (identified spectrum,
peptide, or protein) in separate assays. Each assay is an object of
class
[SingleCellExperiment](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
that allows easy integration with state-of-the-art single-cell analysis
tools. All assays are contained in a single object of class
[QFeatures](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.html).
An overview of the data structure is shown provided in the `scp`
package.

The `scpdata()` function returns a summary table with all currently
available datasets in the package. More information about the data
content and the data collection can be found in the corresponding manual
pages.

## Usage

``` r
scpdata()
```

## Value

A `DataFrame` table containing a summary of the available datasets.

## See also

More information about the data manipulation can be found in the `scp`
package.

## Author

Christophe Vanderaa

## Examples

``` r
## List available datasets and their metadata
scpdata()
#> DataFrame with 39 rows and 15 columns
#>                 title  dataprovider       species taxonomyid      genome
#>           <character>   <character>   <character>  <integer> <character>
#> EH3899  specht2019... SlavovLab ... Homo sapie...       9606          NA
#> EH3900  specht2019... SlavovLab ... Homo sapie...       9606          NA
#> EH3901  dou2019_ly...       MassIVE Homo sapie...       9606          NA
#> EH3902  dou2019_mo...       MassIVE Mus muscul...      10090          NA
#> EH3903  dou2019_bo...       MassIVE Mus muscul...      10090          NA
#> ...               ...           ...           ...        ...         ...
#> EH10674 leduc2025_...        Zenodo Mus muscul...      10090          NA
#> EH10675 leduc2024_...       MassIVE Homo sapie...       9606          NA
#> EH10676 leduc2024_...       MassIVE Homo sapie...       9606          NA
#> EH10677 fulcher202...       MassIVE Homo sapie...       9606          NA
#> EH10678 brunner202...         PRIDE Homo sapie...       9606          NA
#>           description coordinate_1_based    maintainer rdatadateadded
#>           <character>          <integer>   <character>    <character>
#> EH3899  SCP expres...                  1 Christophe...     2020-11-05
#> EH3900  SCP expres...                  1 Christophe...     2020-11-05
#> EH3901  SCP expres...                  1 Christophe...     2020-11-05
#> EH3902  SCP expres...                  1 Christophe...     2020-11-05
#> EH3903  SCP expres...                  1 Christophe...     2020-11-05
#> ...               ...                ...           ...            ...
#> EH10674 Limiting t...                  1 Guillaume ...     2026-11-13
#> EH10675 Single-cel...                  1 Shabnam Za...     2026-11-13
#> EH10676 Single-cel...                  1 Shabnam Za...     2026-11-13
#> EH10677 Single-cel...                  1 Guillaume ...     2026-11-13
#> EH10678 Single-cel...                  1 Léopold Gu...     2026-11-13
#>         preparerclass                                          tags  rdataclass
#>           <character>                                        <AsIs> <character>
#> EH3899        scpdata Experiment...,Expression...,Experiment...,...   QFeatures
#> EH3900        scpdata Experiment...,Expression...,Experiment...,...   QFeatures
#> EH3901        scpdata Experiment...,Expression...,Experiment...,...   QFeatures
#> EH3902        scpdata Experiment...,Expression...,Experiment...,...   QFeatures
#> EH3903        scpdata Experiment...,Expression...,Experiment...,...   QFeatures
#> ...               ...                                           ...         ...
#> EH10674       scpdata      Expression...,MassSpectr...,Proteome,...   QFeatures
#> EH10675       scpdata      Expression...,MassSpectr...,Proteome,...   QFeatures
#> EH10676       scpdata      Expression...,MassSpectr...,Proteome,...   QFeatures
#> EH10677       scpdata      Expression...,MassSpectr...,Proteome,...   QFeatures
#> EH10678       scpdata      Expression...,MassSpectr...,Proteome,...   QFeatures
#>             rdatapath     sourceurl  sourcetype
#>           <character>   <character> <character>
#> EH3899  scpdata/sp... https://sc...         CSV
#> EH3900  scpdata/sp... https://sc...         CSV
#> EH3901  scpdata/do... ftp://mass...    XLS/XLSX
#> EH3902  scpdata/do... ftp://mass...    XLS/XLSX
#> EH3903  scpdata/do... ftp://mass...    XLS/XLSX
#> ...               ...           ...         ...
#> EH10674 records/24... https://ze...         TXT
#> EH10675 records/24... ftp://mass...         TXT
#> EH10676 records/24... ftp://mass...         TXT
#> EH10677 records/24... ftp://mass...         CSV
#> EH10678 records/24... https://ww...         TXT

## Load data using the ExperimentHub interface
hub <- ExperimentHub()

if (FALSE) { # \dontrun{
## Download the data set of interest using ExperimentHub indexing
hub[["EH7711"]]
## Download the same data set using the build-in function
leduc2022()
} # }
```
