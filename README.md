
# Single Cell Proteomics Data Package

This package contains standardized and annotated single-cell proteomics data.

# Installation instruction 

To install the stable version from Bioconductor:

```
if (!requireNamespace("BiocManager"))
	install.packages("BiocManager")
BiocManager::install("scpdata")
```

To install the development version from GitHub:

```
BiocManager::install("UCLouvain-CBIO/scpdata")
```

# Loading a data set 

The available datasets can be accessed through `ExperimentHub`. The 
table with available datasets can be retrieved as follows:

```r
library(scpdata)
scpdata()
```

The desired dataset, for example the `specht2019v2` dataset, can be 
downloaded and accessed using either the bracket selection:

```r
scp <- datasets[["EH3899"]]
```

Or the data loading function:

```r
scp <- specht2019v2()
```

The data description and data collection procedures can be found in 
the corresponding help files that can be accessed using for instance
`?specht2019v2`.

# Available data sets

There are many sources of published singe-cell proteomics data. We 
collected the associated data that can be found under 3 different 
levels:

* **PSM data** is the data obtained after running feature 
  identification and quantification from the raw data. PSM stands for 
  peptide to spectrum match. The rows in the data are PSMs and columns
  are samples.
* **Peptide data** is the expression data where rows are peptides and 
  columns are samples.
* **Protein data** is the expression data where rows are proteins and 
  columns are samples.

Please refer to the documentation and original papers for thorough 
description of the data processing. 
