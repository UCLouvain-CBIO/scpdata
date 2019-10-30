
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

and the data set can be loaded using 

```r
data("DataSetName")
```

For example

```r
data("specht2019")
```

Every data set is documentend. The data collection procedures are described [here](https://github.com/UCLouvain-CBIO/scpdata/inst/README.md), and the data information can be found in the helpfiles. Try

```r
?specht2019
```

## Available data sets

Several articles have been published including single-cell proteomics data. We here distinguish 3 types of data sets: 

* **Raw data** is composed of MS output files and was not processed
* **Peptide data** is the peptide quantitative data obtained after processing the raw data
* **Protein data** is the protein quantitative data obtained after aggregating the peptides belonging to the same protein. 
* $\color{red}{\text{Final expression data after data normalization, imputation, batch correction,... ?}}$

This is an overview of the available data in this package: 

Publ. Date | Data set name | Raw data | Peptide data | Protein data 
---------- | ------------- | -------- | ------------ | ------------
2019-06-09 | specht2019    | Yes      | +/- Yes      | Yes
2019-09-11 | dou2019_1     | Yes      | No           | Yes 
2019-09-11 | dou2019_2     | Yes      | No           | Yes
2019-09-11 | dou2019_3     | Yes      | No           | Yes 

## Data set description 
  
<!-- 
#### Run this and paste output below
desc <- scpdata()$result[, -c(1,2), drop=F]
colnames(desc) <- c("Data", "Description")
knitr::kable(desc) 
-->

|Data       |Description                                                                                                                 |
|:----------|:---------------------------------------------------------------------------------------------------------------------------|
|dou2019_1  |Deep proteome coverage for single cell analysis using nanoPOTS combined with TMT isobaric labeling method (Dou et al. 2019) |
|dou2019_2  |Deep proteome coverage for single cell analysis using nanoPOTS combined with TMT isobaric labeling method (Dou et al. 2019) |
|dou2019_3  |Deep proteome coverage for single cell analysis using nanoPOTS combined with TMT isobaric labeling method (Dou et al. 2019) |
|specht2019 |Quantifying the emergence of macrophage heterogeneity using the SCoPE2 pipeline (Specht et al. 2019)                        |

