
####---- Specht et al. 2019 ---####


# Specht, Harrison, Edward Emmott, Toni Koller, and Nikolai Slavov. 2019. 
# “High-Throughput Single-Cell Proteomics Quantifies the Emergence of Macrophage 
# Heterogeneity.” bioRxiv. https://doi.org/10.1101/665307.

library(MSnbase)
library(tidyr)
setwd("./inst/scripts")

# The data was downloaded from https://scope2.slavovlab.net/docs/data to 
# scpdata/inst/extdata/specht2019


### Experiment metadata

# Create the experiment data for the MSnSet object
expdat <- new("MIAPE",
              title = "High-throughput single-cell proteomics quantifies the emergence of macrophage heterogeneity",
              abstract = "The fate and physiology of individual cells are controlled by networks of proteins. Yet, our ability to quantitatively analyze protein networks in single cells has remained limited. To overcome this barrier, we developed SCoPE2. It integrates concepts from Single-Cell ProtEomics by Mass Spectrometry (SCoPE-MS) with automated and miniaturized sample preparation, substantially lowering cost and hands-on time. SCoPE2 uses data-driven analytics to optimize instrument parameters for sampling more ion copies per protein, thus supporting quantification with improved count statistics. These advances enabled us to analyze the emergence of cellular heterogeneity as homogeneous monocytes differentiated into macrophage-like cells in the absence of polarizing cytokines. We used SCoPE2 to quantify over 2,000 proteins in 356 single monocytes and macrophages in about 85 hours of instrument time, and the quantified proteins allowed us to discern single cells by cell type. Furthermore, the data uncovered a continuous gradient of proteome states for the macrophage-like cells, suggesting that macrophage heterogeneity may emerge even in the absence of polarizing cytokines. Our methodology lays the foundation for quantitative analysis of protein networks at single-cell resolution.",
              url = "https://doi.org/10.1101/665307",
              dateStamp = "2019-07-09",
              name = "Harrison Specht, Edward Emmott, David H. Perlman, Antonius Koller, Nikolai Slavov",
              lab = "Slavov Lab",
              instrumentModel = "Q-Exactive Orbitrap",
              instrumentManufacturer = "Thermo Scientific",
              softwareName = "MaxQuant",
              softwareVersion = "1.6.2.3",
              switchingCriteria = "After a precursor scan from 450 to 1600 m/z at 70,000 resolving power, the top 5 most intense precursor ions with charges 2 to 4 and above the AGC min threshold of 20,000 were isolated for MS2 analysis via a 0.7 Th isolation window",
              ionSource = "ESI",
              ionSourceDetails = "Electrospray voltage was set to 2,200V, applied at the end of the analytical column. To reduce atmospheric background ions and enhance peptide signal to noise ratio, an Active Background Ion Reduction Device (ABIRD, by ESI Source Solutons, LLC, Woburn MA, USA) was used at the nanospray interface. The temperature of ion transfer tube was 250 degrees Celsius and the S-lens RF level set to 80.",
              analyser = "orbitrap",
              analyserDetails = "Precursor ions were accumulated for at most 300ms. Then they were fragmented via HCD at a and the fragments analyzed at 70,000 resolving power. Dynamic exclusion was used with a duration of 30 seconds with a mass tolerance of 10ppm.",
              collisionEnergy = "33 eV (normalized to m/z 500, z=1)")


### Sample metadata

# Get metadata
pdat <- t(read.csv("../extdata/specht2019/Cells.csv", row.names = 1))
# Create the phenotype data for the MSnSet object
pdat <- AnnotatedDataFrame(data = as.data.frame(pdat))


### Peptide data 

# Get data
dat <- read.csv(file = "../extdata/specht2019/Peptides-raw.csv")
rownames(dat) <- dat[, 2]

# Create the expression data for the MSnSet object
edat <- as.matrix(dat[, -c(1,2)])

# Create the feature data for the MSnSet object
fdat <- AnnotatedDataFrame(data = dat[, 1:2])

# Create the MSnSet object
specht2019_peptide <- new("MSnSet", exprs = edat, 
                          featureData = fdat,
                          phenoData = pdat,
                          experimentData = expdat)

# Save data as Rda file
# Note: saving is assumed to occur in "(...)/scpdata/inst/scripts"
save(specht2019_peptide, file = file.path("../../data/specht2019_peptide.rda"),
     compress = "xz", compression_level = 9)



### Protein data

# Get data
dat <- read.csv(file = "../extdata/specht2019/Proteins-processed.csv", row.names = 1)

# Create the expression data for the MSnSet object
edat <- as.matrix(dat[, -ncol(dat)])

# Create the feature data for the MSnSet object
fdat <- AnnotatedDataFrame(data = dat[, ncol(dat), drop = FALSE])

# Create the MSnSet object
specht2019_protein <- new("MSnSet", exprs = edat, 
                          featureData = fdat,
                          phenoData = pdat,
                          experimentData = expdat)

# Save data as Rda file
# Note: saving is assumed to occur in "(...)/scpdata/inst/scripts"
save(specht2019_protein, file = file.path("../../data/specht2019_protein.rda"),
     compress = "xz", compression_level = 9)



####---- Load and process the original peptide data ----####


# MS data
dat0 <- read.table("../extdata/specht2019/ev_updated.txt", header = TRUE, sep = "\t")

# Before formatting data, keep only relevant fields 
intensity.coln <- colnames(dat0)[grepl("^Reporter[.]intensity[.]\\d+$", colnames(dat0))]
.keep <- c("Raw.file", "Modified.sequence", "Sequence", "Length", "Proteins", 
           "Reverse", "Leading.razor.protein", "Potential.contaminant",
           "Gene.names", "Protein.names", "Mass", "Calibrated.retention.time", "PIF", 
           "dart_PEP",  "dart_qval", "Score", intensity.coln)
dat <- dat0[,.keep]


## Correct for isotopic cross contamination

icc <- read.csv("../extdata/specht2019/te269088_lot_correction.csv", row.names = 1)
corrected.ri <- t( solve(icc) %*% t( dat[,intensity.coln] ) )
corrected.ri[corrected.ri < 0.1] <- NA
dat[,intensity.coln] <- corrected.ri


## Filter data

# Remove the reverse hits (from decoy database) and contaminants
dat <- dat[dat$Reverse != "+", ] # sum(dat0$Reverse == "+")/nrow(dat0) * 100 # 20.02% 
dat <- dat[!grepl("^REV", dat$Leading.razor.protein), ] # sum(grepl("^REV", dat0$Leading.razor.protein))/nrow(dat0) * 100 # 20.02% 
dat <- dat[dat$Potential.contaminant != "+", ] # sum(dat0$Potential.contaminant == "+")/nrow(dat0) * 100 # 5.01% 
dat <- dat[dat$PIF > 0.8 & !is.na(dat$PIF), ] # sum(dat0$PIF <= 0.8 | is.na(dat0$PIF))/nrow(dat0) * 100 # 18.21% 

# Remove spectra with poor identification confidence
# This is done with the output of DART-ID
qprots <- unique(dat[dat$dart_qval < 0.01, "Leading.razor.protein"])
dat <- dat[dat$Leading.razor.protein %in% qprots, ]
dat <- dat[dat$dart_PEP < 0.02, ]

# Remove runs with insufficient peptides 
pep.t <- table(dat$Raw.file)
dat <- dat[! dat$Raw.file %in% names(pep.t[pep.t < 300]),]


## Format data

# Change the table to long format by merging the 11 reporter intensities as 2 
# distinct variables (TMT channel and signal intensity)
dat <- as.data.frame(pivot_longer(data = dat, cols = intensity.coln, 
                                  names_to = "Channel", 
                                  values_to = "Reporter.intensity"))
# Deal with duplicate peptides (measured as differently charged ions)
# In case of duplicate peptides in the same run for the same channel, we keep the
# peptide with the lowest PEP
dat <- dat[order(dat$dart_PEP, decreasing = FALSE),]
dat <- dat[!duplicated(dat[,c("Modified.sequence", "Raw.file", "Channel")]), ]
# Create a unique ID for samples 
dat$sample <- paste0(dat$Raw.file, "-", dat$Channel)

## Format the expression data

edat <- pivot_wider(dat, id_cols = "Modified.sequence", names_from = "sample",
                    values_from = "Reporter.intensity", 
                    values_fill = list("Reporter.intensity" = NA))
# 0 intensities are missing values
edat[edat == 0] <- NA
# Add rownames
rown <- edat[,1]
edat <- as.matrix(edat[,-1])
rownames(edat) <- rown$Modified.sequence

# Create the phenotype data 

pdat <- do.call(rbind, lapply(colnames(edat), function(x){
  x <- strsplit(x, "-")[[1]]
  data.frame(run = x[1], channel = x[2])
}))
rownames(pdat) <- colnames(edat)
# Sample metadata 
samp <- read.csv("../extdata/specht2019/annotation_fp60-97.csv", row.names = 1)
batch <- read.csv("../extdata/specht2019/batch_fp60-97.csv", row.names = 1)
ed_specht <- read.csv("../extdata/specht2019/peptides-raw-RI.csv", row.names = 1)
# Add the batch info
pdat <- cbind(pdat, batch[as.character(pdat$run),])
# Add the cell info
pdat <- cbind(pdat, cell_type = do.call(rbind, lapply(rownames(pdat), function(x){
  x <- strsplit(x, "-")[[1]]
  if(grepl("QC", x[1])) return("QC")
  out <- as.character(samp[x[2], paste0("X", x[1])])
  if(length(out) == 0) return(NA)
  return(out)
})))
pdat$channel <- gsub("Reporter[.]intensity[.]", "", pdat$channel)
# Create sample metadata description
mpdat <- data.frame(type = c("characters",
                             "characters",
                             "characters",
                             "characters",
                             "characters",
                             "characters"),
                    labelDescription = c("The ID of the SCoPE-MS run set",
                                         "The TMT channel index within the SCoPE-MS set",
                                         "The chromatographic batch ID",
                                         "The sorting batch ID",
                                         "The digestion batch ID", 
                                         "Phenotype of the sample, determined by experimental design and FACS."))
# Combine all phenotype data in an AnnotatedDataFrame
pdat <- AnnotatedDataFrame(data = pdat, varMetadata = mpdat)

# Create the feature data 

# Keep only variable that are common to all peptides
fdat <- dat[!duplicated(dat$Modified.sequence),
            c("Modified.sequence", "Sequence", "Length", "Proteins", 
              "Gene.names", "Protein.names", "Mass")]
rownames(fdat) <- fdat$Modified.sequence
fdat <- fdat[rownames(edat),]
# Make data descriptioncreate feature metadata description
mfdat <- data.frame(type = c("characters",
                             "characters",
                             "integer",
                             "character",
                             "character",
                             "character",
                             "double"),
                    labelDescription = c("Sequence representation including the post-translational modifications (abbreviation of the modification in brackets before the modified AA). The sequence is always surrounded by underscore characters ('_').",
                                         "The identified amino acid sequence of the peptide.",
                                         "Number of amino acids of the identified peptide.",
                                         "The identifiers of the proteins a particular peptide is associated with.",
                                         "Name(s) of gene(s) the peptide is associated with.",
                                         "Name(s) of the protein(s) the peptide is associated with.",
                                         "The predicted monoisotopic mass of the identified peptide sequence."))
# Combine all feature data in an AnnotatedDataFrame
fdat <- AnnotatedDataFrame(data = fdat, varMetadata = mfdat)

# Create the experimental metadata
# Collected from article
expdat <- new("MIAPE",
              title = "High-throughput single-cell proteomics quantifies the emergence of macrophage heterogeneity",
              abstract = "The fate and physiology of individual cells are controlled by networks of proteins. Yet, our ability to quantitatively analyze protein networks in single cells has remained limited. To overcome this barrier, we developed SCoPE2. It integrates concepts from Single-Cell ProtEomics by Mass Spectrometry (SCoPE-MS) with automated and miniaturized sample preparation, substantially lowering cost and hands-on time. SCoPE2 uses data-driven analytics to optimize instrument parameters for sampling more ion copies per protein, thus supporting quantification with improved count statistics. These advances enabled us to analyze the emergence of cellular heterogeneity as homogeneous monocytes differentiated into macrophage-like cells in the absence of polarizing cytokines. We used SCoPE2 to quantify over 2,000 proteins in 356 single monocytes and macrophages in about 85 hours of instrument time, and the quantified proteins allowed us to discern single cells by cell type. Furthermore, the data uncovered a continuous gradient of proteome states for the macrophage-like cells, suggesting that macrophage heterogeneity may emerge even in the absence of polarizing cytokines. Our methodology lays the foundation for quantitative analysis of protein networks at single-cell resolution.",
              url = "https://doi.org/10.1101/665307",
              dateStamp = "2019-07-09",
              name = "Harrison Specht, Edward Emmott, David H. Perlman, Antonius Koller, Nikolai Slavov",
              lab = "Slavov Lab",
              instrumentModel = "Q-Exactive Orbitrap",
              instrumentManufacturer = "Thermo Scientific",
              softwareName = "MaxQuant",
              softwareVersion = "1.6.2.3",
              switchingCriteria = "After a precursor scan from 450 to 1600 m/z at 70,000 resolving power, the top 5 most intense precursor ions with charges 2 to 4 and above the AGC min threshold of 20,000 were isolated for MS2 analysis via a 0.7 Th isolation window",
              ionSource = "ESI",
              ionSourceDetails = "Electrospray voltage was set to 2,200V, applied at the end of the analytical column. To reduce atmospheric background ions and enhance peptide signal to noise ratio, an Active Background Ion Reduction Device (ABIRD, by ESI Source Solutons, LLC, Woburn MA, USA) was used at the nanospray interface. The temperature of ion transfer tube was 250 degrees Celsius and the S-lens RF level set to 80.",
              analyser = "orbitrap",
              analyserDetails = "Precursor ions were accumulated for at most 300ms. Then they were fragmented via HCD at a and the fragments analyzed at 70,000 resolving power. Dynamic exclusion was used with a duration of 30 seconds with a mass tolerance of 10ppm.",
              collisionEnergy = "33 eV (normalized to m/z 500, z=1)")

# Create the MSnSet object
specht2019_peptide2 <- new("MSnSet", exprs = edat, 
                           featureData = fdat,
                           phenoData = pdat,
                           experimentData = expdat)
  
# Save data as Rda file
# Note: saving is assumed to occur in "(...)/scpdata/inst/scripts"
save(specht2019_peptide2, file = file.path("../../data/specht2019_peptide2.rda"),
     compress = "xz", compression_level = 9)







