
# Single Cell Proteomics Data Package

This package contains formatted and annotated single cell proteomics data.

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
scpdata("DataSetName", type = "dataType")
```

For example

```r
sc <- scpdata("specht2019", type = "peptide")
```

The data description and data collection procedures can be found in the corresponding helpfiles. Try

```r
?specht2019
```

# Available data sets

Several articles have been published including single-cell proteomics data. We here distinguish 3 types of data sets: 

* **Raw data** is composed of MS output files and was not processed
* **Peptide data** is the peptide quantitative data obtained after processing the raw data
* **Protein data** is the protein quantitative data obtained after aggregating the peptides belonging to the same protein. 

The amount of processing of each type of data might differ from data set to data set. Please refer to the documentation and original papers for thorough description of the data processing. 

This is an overview of the available data in this package: 

Publ. Date | Data set name | Raw data | Peptide data | Protein data 
---------- | ------------- | -------- | ------------ | ------------
2019-06-09 | specht2019    | No       | Yes          | Yes
2019-09-11 | dou2019_hela     | No       | No           | Yes 
2019-09-11 | dou2019_boosting     | No       | No           | Yes
2019-09-11 | dou2019_mouse     | No       | No           | Yes 
2018-08-24 | specht2018    | No       | Yes          | No

# Data set description 
  
<!-- 
#### Run this and paste output below
desc <- scpdata()$result[, -c(1,2), drop=F]
colnames(desc) <- c("Data", "Description")
knitr::kable(desc) 
-->

Data               |Description                                                                                                                                                                        |
|:------------------|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|dou2019_hela_protein  |High-Throughput Single Cell Proteomics Enabled by Multiplex Isobaric Labeling in a Nanodroplet Sample Preparation Platform: HeLa digests (Dou et al. 2019)                         |
|dou2019_boosting_protein  |High-Throughput Single Cell Proteomics Enabled by Multiplex Isobaric Labeling in a Nanodroplet Sample Preparation Platform: testing boosting ratios (Dou et al. 2019)              |
|dou2019_mouse_protein  |High-Throughput Single Cell Proteomics Enabled by Multiplex Isobaric Labeling in a Nanodroplet Sample Preparation Platform: profiling of murine cell populations (Dou et al. 2019) |
|specht2018_peptide |mPOP SCoPE-MS Master Mix 20180824 (Specht et al. 2018)                                                                                                                             |
|specht2019_peptide |Quantifying the emergence of macrophage heterogeneity using the SCoPE2 pipeline (Specht et al. 2019)                                                                               |
|specht2019_protein |Quantifying the emergence of macrophage heterogeneity using the SCoPE2 pipeline (Specht et al. 2019)                                        