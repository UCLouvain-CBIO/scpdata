
# Single Cell Proteomics Data Package

This package contains standardized and annotated single-cell proteomics data.

# Installation instruction 

```r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("UCLouvain-CBIO/scpdata")
```

# Loading a data set 

The available datasets can be accessed through the `ExperimentHub` API. The list
of available datasets can be retrieved as follows.

```r
library(ExperimentHub)
eh <- ExperimentHub()
(datasets <- query(eh, "scpdata"))
```

The desired dataset, for example the `specht2019v2` dataset, can be downloaded 
and accessed using either the bracket selection:

```r
scp <- datasets[["EH????"]]
```

Or the data loading function:

```r
scp <- specht2019v2(metadata = FALSE)
```

The data description and data collection procedures can be found in the 
corresponding helpfiles that can be accessed using for instance `?specht2019v2`.

# Available data sets

There are many sources of published singe-cell proteomics work. We collected the
associated data that can be found under 3 different levels:

* **PSM data** is the data obtained after running feature identification and 
quantification from the raw data. Rows are spectral features and columns are 
samples.
* **Peptide data** is the expression data where rows are peptides and columns 
are samples
* **Protein data** is the expression data where rows are proteins and columns 
are samples

The amount of processing of each type of data might differ from data set to data 
set. Please refer to the documentation and original papers for thorough 
description of the data processing. 

Meta-information about the available datasets in `scpdata` can be fetched with 
the following command: 

```r
mcols(datasets)
```