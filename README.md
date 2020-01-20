
# Single Cell Proteomics Data Package

This package contains standardized and annotated single-cell, or close to single-cell, proteomics data.

# Installation instruction 

```r
if (!requireNamespace("devtools", quietly=TRUE))
    install.packages("devtools")
devtools::install_github("UCLouvain-CBIO/scpdata")
```

# Loading a data set 

The available data sets can be listed using 

```r
scpdata()
```

and the data set can be accessed using 

```r
data("DataSetName")
# Example
data("specht2018_peptide")
```

The data description and data collection procedures can be found in the corresponding helpfiles. Try

```r
?specht2018_peptide
```

# Available data sets

Several articles have been published including single-cell proteomics data. We here distinguish 3 types of data sets: 

* **Raw data** is the unprocessed mass spectrometry data
* **Peptide data** is the expression data where rows are peptides and columns are samples
* **Protein data** is the expression data where rows are proteins and columns are samples

The amount of processing of each type of data might differ from data set to data set. Please refer to the documentation and original papers for thorough description of the data processing. 

This is an overview of the available data in this package: 

Publ. Date | Data set name    | Data type  | Raw data                                                                   | Peptide data | Protein data 
---------- | ---------------- | ---------- | -------------------------------------------------------------------------- | ------------ | -------
2019-09-11 | dou2019_hela     | multiplex  | [ftp link](ftp://massive.ucsd.edu/MSV000084110)                            | No           | Yes 
2019-09-11 | dou2019_boosting | multiplex  | [ftp link](ftp://massive.ucsd.edu/MSV000084110)                            | No           | Yes
2019-09-11 | dou2019_mouse    | multiplex  | [ftp link](ftp://massive.ucsd.edu/MSV000084110)                            | No           | Yes 
2019-06-09 | specht2019v1     | multiplex  | [ftp link](ftp://massive.ucsd.edu/MSV000083945)                            | Yes          | Yes
2018-08-24 | specht2018       | multiplex  | [ftp link](ftp://massive.ucsd.edu/MSV000082841)                            | Yes          | No
2018-02-28 | zhu2018NC_hela   | label free | [ftp link](ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2018/01/PXD006847) | Yes          | No
2018-02-28 | zhu2018NC_islets | label free | [ftp link](ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2018/01/PXD006847) | Yes          | No
2018-01-10 | zhu2018MCP       | label free | [ftp link](ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2018/07/PXD008844) | Yes          | No

# Data set description 
  
<!-- 
#### Run this and paste output below
desc <- scpdata()$result[, -c(1,2), drop=F]
colnames(desc) <- c("Data", "Description")
knitr::kable(desc) 
-->

|Data                     |Description                                                                   |
|:------------------------|:-----------------------------------------------------------------------------|
|dou2019_boosting_protein |FACS + nanoPOTS + TMT: testing boosting ratios (Dou et al. 2019)              |
|dou2019_hela_protein     |FACS + nanoPOTS + TMT: HeLa digests (Dou et al. 2019)                         |
|dou2019_mouse_protein    |FACS + nanoPOTS + TMT: profiling of murine cell populations (Dou et al. 2019) |
|specht2018_peptide       |SCoPE-MS + mPOP lysis upgrade: Master Mix 20180824 (Specht et al. 2018)       |
|specht2019v1_peptide     |FACS + SCoPE2: macrophages vs monocytes (Specht et al. 2019)                  |
|specht2019v1_protein     |FACS + SCoPE2: macrophages vs monocytes (Specht et al. 2019)                  |
|zhu2018MCP_peptide       |laser dissection + nanoPOTS: (Zhu et al. 2018, Mol Cell Proteomics)           |
|zhu2018NC_hela_peptide   |nanoPOTS: HeLa dilutions (Zhu et al. 2018, Nat Comm)                          |
|zhu2018NC_islets_peptide |laser dissection + nanoPOTS: T1D islets (Zhu et al. 2018, Nat Comm)           |



