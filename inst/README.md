
# Methodology used for collecting the data


1. Identify the data source and the annotations from the article, ask authors 
   for data if needed. 
2. Create a new R script that converts the data using the `scp` framework
3. Add data documentation and the data collection procedure in 
   `scpdata/R/data.R`
4. Add the dataset information in the `make-metadata.R` and update the 
   `metadata.csv`
5. Test the new data set by running 
  `ExperimentHubData::makeExperimentHubMetadata("scpdata")`. 
6. Contact Bioc team and upload Rda to Microsoft Azure: 
   hubs@bioconductor.org
   Then upload the files using the provided link:
   `azcopy copy --recursive scpdata/ '<sas-url>'`
   or more easily, use the R script. See the 
   [help page](https://bioconductor.org/packages/devel/bioc/vignettes/HubPub/inst/doc/CreateAHubPackage.html#uploading-data-to-microsoft-azure-genomic-data-lake)
   for more information
7. Compile the documentation with roxygen2 and check package.
8. Update the NEWS.md file and bump package version

