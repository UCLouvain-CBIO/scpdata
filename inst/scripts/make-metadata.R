

####---- Make the ExperimentHub metadata table ----####

meta <- data.frame(
  Title = c(paste0("FACS + SCoPE2: macrophages vs monocytes (Specht et al. ",
                   "2019 - version 2)")),
  Description = c(paste0("Single cell proteomics data acquired by the Slavov ",
                         "Lab. This is the version 2 of the data released in ", 
                         "December 2019. It contains quantitative information ", 
                         "of macrophages and monocytes at PSM, peptide and ", 
                         "protein level.")),
  BiocVersion = "3.11",
  Genome = NA_character_,
  SourceType = c("CSV"),
  SourceUrl = c(paste0("https://drive.google.com/drive/folders/1VzBfmNxziRYqayx3SP-cOe2gu129Obgx")),
  SourceVersion=NA_character_,
  Species=NA_character_,
  TaxonomyId=NA_integer_,
  Coordinate_1_based=TRUE,
  DataProvider=NA_character_,
  aintainer=NA_character_,
  RDataClass=NA_character_,
  DispatchClass=NA_character_,
  RDataPath=NA_character_,
  Tags=NA_character_,
  Notes=NA_character_)


