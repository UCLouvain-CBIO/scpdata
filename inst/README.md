

# Methodology used for collecting the data


1. Identify the data source and the annotations from the article, ask authors for data if needed
2. Create a new R script that converts the data to `MSnSet` objects
3. Add data documentation in `scpdata/R/data.R`
4. Update the `scpdata/README.md` file


# Dou et al. 2019

*Dou, Maowei, Geremy Clair, Chia-Feng Tsai, Kerui Xu, William B. Chrisler, Ryan L. Sontag, Rui Zhao, et al. 2019. “High-Throughput Single Cell Proteomics Enabled by Multiplex Isobaric Labeling in a Nanodroplet Sample Preparation Platform.” Analytical Chemistry, September. https://doi.org/10.1021/acs.analchem.9b03349.*

The article contains 3 SCP data sets available on the [ACS Publications](https://pubs.acs.org/doi/10.1021/acs.analchem.9b03349) website:

* Supplementary data set 1, raw data for HeLa digest (`.xlsx`)
* Supplementary data set 2, raw data for testing boosting ratios (`.xlsx`)
* Supplementary data set 3, raw data for isobaric labelling-based single cell quantification and bulk-scale label free quantification (`.xlsx`)

## Supplementary data set 1

We arbitrarily call this data set `dou2019_1`. The `.xlsx` spreadsheet contrains the following sheets:

* `0 - Description`
* `1 - Run 1 raw data`: This table contains the raw data for the Run 1 only. It correspond to the assembly of MSGF+ identifications and MASIC reporter assemblies. The data was then isotope corrected and sum rolled-up to the protein level. Contaminant and reverse hit were removed from this table. Median values are indicated at the bottom of the table
* `2 - Run 1 processed data`: This table contains the raw data for the Run 1 alone. The protein data was log2 transformed, median normalized, the "sva::ComBat" function was used to remove TMT-set-dependent Batch effects. 
* `3 - Run 2 raw data`:	This table contains the raw data for the Run 2 only. It correspond to the assembly of MSGF+ identifications and MASIC reporter assemblies. The data was then isotope corrected and sum rolled-up to the protein level. Contaminant and reverse hit were removed from this table. Median values are indicated at the bottom of the table
* `4 - Run 2 processed data`:	This table contains the raw data for the Run 2 alone. The protein data was log2 transformed, median normalized, the "sva::ComBat" function was used to remove TMT-set-dependent Batch effects. 
* `5 - Run 1 and 2 raw data`:	This table contains the raw data for the *[Run 1 only]* (the authors probably mean Run 1 and Run 2). It correspond to the assembly of MSGF+ identifications and MASIC reporter assemblies. The data was then isotope corrected and sum rolled-up to the protein level. Contaminant and reverse hit were removed from this table. Median values are indicated at the bottom of the table
* `6 - Run 1 and 2 processed data`:	This table contains the raw data for the Run 1 and 2 together. The protein data was log2 transformed, median normalized, the "sva::ComBat" function was used to remove TMT-set-dependent Batch effects. The crossdataset missingness is indicated at the right of the table

We only used the sheet 5 (`Run 1 and 2 raw data`) containing the combined "raw" data for the two runs. This data set is converted to an `MSnSet` object in the script `scpdata/inst/script/dou2019.R`.

## Supplementary data set 2

We arbitrarily call this data set `dou2019_2`. The `.xlsx` spreadsheet contrains the following sheets:

* `0 - Description`
* `01 - No Boost raw data`:	This table contains the raw data for the Run 1 and 2 of the no boost. It correspond to the assembly of MSGF+ identifications and MASIC reporter assemblies. The data was then isotope corrected and sum rolled-up to the protein level. Contaminant and reverse hit were removed from this table.
* `02 - No Boost processed data`:	This table contains the raw data for the Run 1 and 2 of the no boost. The protein data was log2 transformed, median normalized, the "sva::ComBat" function was used to remove TMT-set-dependent Batch effects. 
* `03 - 5ng boost raw data`:	This table contains the raw data for the Run 1 and 2 of the 5ng boost. It correspond to the assembly of MSGF+ identifications and MASIC reporter assemblies. The data was then isotope corrected and sum rolled-up to the protein level. Contaminant and reverse hit were removed from this table.
* `04 - 5ng boost processed data`:	This table contains the raw data for the Run 1 and 2 of the 5ng boost. The protein data was log2 transformed, median normalized, the "sva::ComBat" function was used to remove TMT-set-dependent Batch effects. 
* `05 - 50ng boost raw data`:	This table contains the raw data for the Run 1 and 2 of the 50ng boost. It correspond to the assembly of MSGF+ identifications and MASIC reporter assemblies. The data was then isotope corrected and sum rolled-up to the protein level. Contaminant and reverse hit were removed from this table.
* `06 - 50ng processed data`:	This table contains the raw data for the Run 1 and 2 of the 50ng boost. The protein data was log2 transformed, median normalized, the "sva::ComBat" function was used to remove TMT-set-dependent Batch effects. 

We combined the data from sheet 01, 03, and 05, that is all tables containing the "raw" data. This corresponds to merging the data sets with no boosting, 5ng boosting, and 50 ng boosting. The boosting quantities are kept as feature data. Combining data is performed by matching the proteins beetween boosting batches. The data set is formated to an `MSnSet` object in the script `scpdata/inst/script/dou2019.R`.

## Supplementary data set 3

We arbitrarily call this data set `dou2019_3`. The `.xlsx` spreadsheet contrains the following sheets:

* `00 - Description`
* `01 - Raw sc protein data` Contains the assembly of MSGF+ identifications and MASIC assemblies isotope corrected and sum rolled-up to the protein level. Contaminant and reverse hit were removed from this table
* `02 - Processed sc protein data` Only the proteins with at least 2 unique peptides were conserved for further analysis the protein data was log2 transformed, filtered for outliers using the pmartR method (PMID: 30638385 ), median normalized, only proteins with at least 60% of values within a cell type were conserved for further analysis, the "sva::ComBat" function was used to remove TMT-set-dependent Batch effects. Paired two-tailed heteroscedastic T-tests were performed in order to establish protein enriched in a given condition and Z-scores were calculated for heatmap visualization.
* `03 - Bulk Proteomics` The bulk proteomics was generated on a lysate of each cell type in 5 technical replicates. iBAQ intensities were used the data was log transformed, median normalized, imputed using a normal distribution of width 0.5 at 1.8 standard deviation away from the median of the data distribution (as described in PMID: 27348712).
* `04 - C10_vs_Raw` This sheet contain the benchmaring of the single cell results to the bulk data for the C10 vs Raw comparison.
* `05 - SVEC_vs_Raw` This sheet contain the benchmaring of the single cell results to the bulk data for the SVEC vs Raw comparison.
* `06 - C10_vs_SVEC` This sheet contain the benchmaring of the single cell results to the bulk data for the C10 vs SVEC comparison.

We parsed the sheet 02 (`Raw sc protein data`) to an MSnSet without further modification (see script `scpdata/inst/script/dou2019.R`). 


# Specht et al. 2018

*Specht, Harrison, Guillaume Harmange, David H. Perlman, Edward Emmott, Zachary Niziolek, Bogdan Budnik, and Nikolai Slavov. 2018. “Automated Sample Preparation for High-Throughput Single-Cell Proteomics.” bioRxiv. https://doi.org/10.1101/399774.*

The MaxQuant output files are available on the [Google Drive](https://drive.google.com/drive/folders/19YG70I52DH5yETcZagdUjNZWNPs0JXVr) account of the SlavovLab (see also [here](http://slavovlab.net/mPOP/index.html). We used the following file:

* `evidence.txt`: the raw output of the MaxQuant software. Documentation of the fields can be found on the [MaxQuant website](http://www.coxdocs.org/doku.php?id=maxquant:table:evidencetable)

The `evidence.txt` file is loaded from the Google Drive folder to R and is processed as follows: 

* Contaminants and reverse hits are removed
* Peptides with PIF (parent ion fraction) smaller than 0.8 and PEP (posterior error probability) greater then 0.02 are removed
* Experiment sets with less than 300 identified peptides are discarded. Note that all experiment sets contained more than 300 identified peptides.
* Some peptides were identified twice within the same experiment set because of different charge states. The peptide identification with lowest PEP is kept and the other(s) discarded. 
* Intensity data is formated to a feature (peptide) x sample (combination of experiment run and TMT channel) matrix. 
* Peptide information is gathered in a feature data frame. PEP, PIF, Score, and retention time are aggregated by taking the median. When the mass is differing among the same peptide sequence, only the mass for the unmodified peptide is kept.
* Sample information is gathered in a phenotype data frame.
* All information is stored in an `MSnSet` object

These processing steps can be found in `scpdata/inst/script/specht2018.R`

