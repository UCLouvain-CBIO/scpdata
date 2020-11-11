
# Methodology used for collecting the data


1. Identify the data source and the annotations from the article, ask authors 
for data if needed. 
2. Create a new R script that converts the data using the `scp` framework
3. Add data documentation and the data collection procedure in 
`scpdata/R/data.R`
4. Add the dataset information in the `make-metadata.R` and update the 
`metadata.csv`
5. Test the new data set. Once ok, upload it to AWS S3 bucket. 
