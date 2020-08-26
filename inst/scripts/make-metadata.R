

####---- Make the ExperimentHub metadata table ----####

meta <- data.frame(
  Title = c(paste0("FACS + SCoPE2: macrophages vs monocytes (Specht et al. ",
                   "2019 - version 2)")),
  Description = c(paste0("Single cell proteomics data acquired by the Slavov ",
                         "Lab. This is the version 2 of the data released in ", 
                         "December 2019. It contains quantitative information ", 
                         "of monocytes (U-937) and macrophages at PSM, peptide ",
                         "and protein level.")),
  BiocVersion = "3.11",
  Genome = NA_character_,
  SourceType = c("CSV"),
  SourceUrl = c("https://scope2.slavovlab.net/docs/data"),
  SourceVersion = c(NA_character_),
  Species = "Homo sapiens",
  TaxonomyId = 9606,
  Coordinate_1_based = TRUE,
  DataProvider = "SlavovLab",
  Maintainer = "Christophe Vanderaa <christophe.vanderaa@uclouvain.be>",
  RDataClass = "QFeatures",
  DispatchClass = "Rda",
  RDataPath = "scpdata/specht2019v2.Rda",
  Notes=NA_character_)

write.csv(meta, "inst/extdata/metadata.csv")

