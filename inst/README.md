
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
6. Contact Bioc team and upload Rda to AWS S3 bucket. To do so, 
   connect to amazon:
7. Compile the documentation with roxygenize.
8. Update the NEWS.md file and bump package version


```
## Login
aws configure --profile AnnotationContributor
## Upload data
aws --profile AnnotationContributor s3 cp scpdata s3://annotation-contributor/scpdata --recursive --acl public-read
## Check
aws --profile AnnotationContributor s3 ls s3://annotation-contributor/scpdata
```
