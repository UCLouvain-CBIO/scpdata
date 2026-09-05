# Single-Cell Proteomics Data Package

The `scpdata` package distributes mass spectrometry-based single-cell
proteomics datasets. The datasets were collected from published work and
formatted to a standardized data framework. The `scp` frameworks stores
the expression data for different MS levels (identified spectrum,
peptide, or protein) in separate assays. Each assay is an object of
class
[SingleCellExperiment::SingleCellExperiment](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
that allows easy integration with state-of-the-art single-cell analysis
tools. All assays are contained in a single object of class
[QFeatures::QFeatures](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.html).
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
#> DataFrame with 32 rows and 15 columns
#>                 title  dataprovider       species taxonomyid      genome
#>           <character>   <character>   <character>  <integer> <character>
#> EH3899  specht2019... SlavovLab ... Homo sapie...       9606          NA
#> EH3900  specht2019... SlavovLab ... Homo sapie...       9606          NA
#> EH3901  dou2019_ly...       MassIVE Homo sapie...       9606          NA
#> EH3902  dou2019_mo...       MassIVE Mus muscul...      10090          NA
#> EH3903  dou2019_bo...       MassIVE Mus muscul...      10090          NA
#> ...               ...           ...           ...        ...         ...
#> EH9610  hu2023_K56...       MassIVE Homo sapie...       9606          NA
#> EH9611  hu2023_ooc...       MassIVE Homo sapie...       9606          NA
#> EH9627        ai2025a       MassIVE Homo sapie...       9606          NA
#> EH10516     leduc2025        Zenodo Mus muscul...      10090          NA
#> EH10517        wu2026         PRIDE Homo sapie...       9606          NA
#>           description coordinate_1_based    maintainer rdatadateadded
#>           <character>          <integer>   <character>    <character>
#> EH3899  SCP expres...                  1 Christophe...     2020-11-05
#> EH3900  SCP expres...                  1 Christophe...     2020-11-05
#> EH3901  SCP expres...                  1 Christophe...     2020-11-05
#> EH3902  SCP expres...                  1 Christophe...     2020-11-05
#> EH3903  SCP expres...                  1 Christophe...     2020-11-05
#> ...               ...                ...           ...            ...
#> EH9610  Single-cel...                  1 Enes Sefa ...     2024-11-13
#> EH9611  Single-cel...                  1 Enes Sefa ...     2024-11-13
#> EH9627  Single-cel...                  1 Laurent Ga...     2025-03-06
#> EH10516 Principles...                  1 Enes Sefa ...     2026-08-18
#> EH10517 Label-free...                  1 Léopold Gu...     2026-09-04
#>         preparerclass                                          tags
#>           <character>                                        <AsIs>
#> EH3899        scpdata Experiment...,Expression...,Experiment...,...
#> EH3900        scpdata Experiment...,Expression...,Experiment...,...
#> EH3901        scpdata Experiment...,Expression...,Experiment...,...
#> EH3902        scpdata Experiment...,Expression...,Experiment...,...
#> EH3903        scpdata Experiment...,Expression...,Experiment...,...
#> ...               ...                                           ...
#> EH9610        scpdata      Expression...,MassSpectr...,Proteome,...
#> EH9611        scpdata      Expression...,MassSpectr...,Proteome,...
#> EH9627        scpdata      Expression...,MassSpectr...,Proteome,...
#> EH10516       scpdata      Expression...,MassSpectr...,Proteome,...
#> EH10517       scpdata      Expression...,MassSpectr...,Proteome,...
#>            rdataclass     rdatapath     sourceurl  sourcetype
#>           <character>   <character>   <character> <character>
#> EH3899      QFeatures scpdata/sp... https://sc...         CSV
#> EH3900      QFeatures scpdata/sp... https://sc...         CSV
#> EH3901      QFeatures scpdata/do... ftp://mass...    XLS/XLSX
#> EH3902      QFeatures scpdata/do... ftp://mass...    XLS/XLSX
#> EH3903      QFeatures scpdata/do... ftp://mass...    XLS/XLSX
#> ...               ...           ...           ...         ...
#> EH9610  SingleCell... scpdata/hu... ftp://mass...         TXT
#> EH9611  SingleCell... scpdata/hu... ftp://mass...         TXT
#> EH9627      QFeatures scpdata/ai... ftp://mass...         TXT
#> EH10516     QFeatures records/21... https://ze...         CSV
#> EH10517     QFeatures records/22... https://ww...         TXT

## Load data using the ExperimentHub interface
hub <- ExperimentHub()

if (FALSE) { # \dontrun{
## Download the data set of interest using ExperimentHub indexing
hub[["EH7711"]]
## Download the same data set using the build-in function
leduc2022()
} # }
```
