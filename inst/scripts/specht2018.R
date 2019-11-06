
####---- Specht et al. 2018 ---####


# Specht, Harrison, Guillaume Harmange, David H. Perlman, Edward Emmott, Zachary 
# Niziolek, Bogdan Budnik, and Nikolai Slavov. 2018. “Automated Sample 
# Preparation for High-Throughput Single-Cell Proteomics.” bioRxiv. 
# https://doi.org/10.1101/399774.

library(MSnbase)
library(tidyr)

# Download data
# MS data
dat0 <- read.table("https://drive.google.com/uc?export=download&id=19o-vtyKOmVmlOoiPH4kri3mB4N3zUUcC", 
                   header = TRUE, sep = "\t")
dat <- dat0
intensity.coln <- colnames(dat)[grepl("^Reporter[.]intensity[.]\\d+$", colnames(dat))]

# Sample data 
# TODO ask for exp design

# Data QC

# Convert the original data to wide format
datw <- pivot_longer(data = dat0, cols = intensity.coln, 
                    names_to = "Channel", 
                    values_to = "Reporter.intensity")
# Plot the RI per channel to confirm empty and carrier channels
axis.txt <- array(data = paste0("TMT", 1:11), dim = 11,
                  dimnames = list(intensity.coln))
ggplot(data = datw, aes(x = Channel, y = Reporter.intensity)) +
  geom_violin() + scale_y_log10() +
  scale_x_discrete(labels = axis.txt,
                   limits = paste0("Reporter.intensity.", 0:10))
# CCL: Channel 1 and 2 are carrier and channel 11 is empty
# Plot the RI per run to identify failed runs
ggplot(data = datw, aes(x = Raw.file, y = Reporter.intensity)) +
  geom_violin() + scale_y_log10() + 
  theme(axis.text.x = element_text(angle=90))
# CCL: no failed run


# Correct for isotopic cross correlation
# TODO ask for cross contamination table

## Filter data

# Remove the reverse hits (from decoy database) and contaminants
dat <- dat[dat$Reverse != "+", ] # sum(dat0$Reverse == "+")/nrow(dat0) * 100 # 20.02% 
dat <- dat[!grepl("^REV", dat$Leading.razor.protein), ] # sum(grepl("^REV", dat0$Leading.razor.protein))/nrow(dat0) * 100 # 20.02% 
dat <- dat[dat$Potential.contaminant != "+", ] # sum(dat0$Potential.contaminant == "+")/nrow(dat0) * 100 # 5.01% 
dat <- dat[dat$PIF > 0.8 & !is.na(dat$PIF), ] # sum(dat0$PIF <= 0.8 | is.na(dat0$PIF))/nrow(dat0) * 100 # 18.21% 

# Remove spectra with poor identification confidence
# TODO missing FDR ...
dat <- dat[dat$PEP < 0.02, ] # sum(dat0$PEP > 0.02)/nrow(dat0) * 100 # 51.00 % 

# Remove runs with insufficient peptides 
pep.t <- table(dat$Raw.file)
dat <- dat[! dat$Raw.file %in% names(pep.t[pep.t < 300]),]

# Other filters based on carrier ratio
# TODO when experiment design data is available

## Format data

# Change the table to long format by merging the 11 reporter intensities as 2 
# distinct variables (TMT channel and signal intensity)
dat <- pivot_longer(data = dat, cols = intensity.coln, 
                    names_to = "Channel", 
                    values_to = "Reporter.intensity")

# Deal with duplicate peptides (measured as differently charged ions)
# In case of duplicate peptides in the same run for the same channel, we keep the
# peptide with the lowest PEP
dat <- dat[order(dat$PEP, decreasing = TRUE),]
dat <- dat[!duplicated(dat[,c("Sequence", "Raw.file", "Channel")]), ]

# Create a unique ID for samples 
dat$sample <- paste0(dat$Raw.file, "-", dat$Channel)

# Format the expression data
edat <- pivot_wider(dat, id_cols = "Sequence", names_from = "sample",
                    values_from = "Reporter.intensity", 
                    values_fill = list("Reporter.intensity" = NA))
rown <- edat[,1]
edat <- as.matrix(edat[,-1])
rownames(edat) <- rown$Sequence

# Create the phenotype data 
pdat <- do.call(rbind, lapply(colnames(edat), function(x){
  x <- strsplit(x, "-")[[1]]
  data.frame(run = x[1], channel = x[2])
}))
rownames(pdat) <- colnames(edat)
# Create sample metadata description
mpdat <- data.frame(type = c("characters",
                             "characters"),
                    labelDescription = c("the id of the SCoPE-MS run set",
                                         "TMT channel id within the SCoPE-MS set"))
# Combine all phenotype data in an AnnotatedDataFrame
pdat <- AnnotatedDataFrame(data = pdat, varMetadata = mpdat)

# Create the feature data 
fdat <- do.call(rbind, lapply(rownames(edat), function(x){
  .idx <- dat$Sequence == x
  .sub <- dat[.idx, , drop = FALSE]
  # Add common field 
  line <- unique(.sub[, c("Sequence", "Length", "Proteins", 
                          "Gene.names", "Protein.names"),])
  # Add peptide count, that is in how many samples was the peptide found
  line$Count <- sum(.idx)
  # Take the smallest theoretical mass (mass for the unmodifiied peptide)
  line$Mass <- min(.sub$Mass)
  # Average the remaining fields fields
  line <- data.frame(line, 
                     as.list(apply(.sub[,c("Calibrated.retention.time", "PIF", 
                                           "PEP", "Score")], 2, median, na.rm = TRUE)))
  return(line)
}))
rownames(fdat) <- rownames(edat)
# Make data descriptioncreate feature metadata description
mfdat <- data.frame(type = c("characters",
                             "integer",
                             "character",
                             "character",
                             "character",
                             "integer",
                             "double",
                             "double",
                             "double",
                             "double",
                             "double"),
                    labelDescription = c("The identified amino acid sequence of the peptide.",
                                         "Number of amino acids of the identified peptide.",
                                         "The identifiers of the proteins a particular peptide is associated with.",
                                         "Name(s) of gene(s) the peptide is associated with.",
                                         "Name(s) of the protein(s) the peptide is associated with.",
                                         "The number of samples where the peptide was identified.",
                                         "The predicted monoisotopic mass of the identified peptide sequence.",
                                         "The recalibrated retention time in minutes in the elution profile of the precursor ion. The value is the median over all samples where the peptide was identified.",
                                         "Short for Parent Ion Fraction; indicates the fraction the target peak makes up of the total intensity in the inclusion window. The value is the median over all samples where the peptide was identified.",
                                         "Posterior Error Probability of the identification. This value essentially operates as a p-value, where smaller is more significant. The value is the median over all samples where the peptide was identified.",
                                         "Andromeda score for the best associated MS/MS spectrum. The value is the median over all samples where the peptide was identified."))
# Combine all feature data in an AnnotatedDataFrame
fdat <- AnnotatedDataFrame(data = fdat, varMetadata = mfdat)

# Create the experimental metadata
expdat <- new("MIAPE",
              title = "Automated sample preparation for high-throughputsingle-cell proteomics",
              abstract = "A major limitation to applying quantitative LC-MS/MS proteomics to small samples, suchas single cells, are the losses incured during sample cleanup. To relieve this limitation, we de-veloped a Minimal ProteOmic sample Preparation (mPOP) method for culture-grown mam-malian cells. mPOP obviates cleanup and thus eliminates cleanup-related losses while expe-diting sample preparation and simplifying its automation.  Bulk SILAC samples processedby mPOP or by conventional urea-based methods indicated that mPOP results in completecell lysis and accurate relative quantification. We integrated mPOP lysis with the Single CellProtEomics by Mass Spectrometry (SCoPE-MS) sample preparation, and benchmarked thequantification  of  such  samples  on  a  Q-exactive  instrument.   The  results  demonstrate  lownoise and high technical reproducibility. Then, we FACS sorted single U-937, HEK-293, andmouse ES cells into 96-well plates and analyzed them by automated mPOP and SCoPE-MS.The quantified proteins enabled separating the single cells by cell-type and cell-division-cyclephase.",
              url = "http://dx.doi.org/10.1101/399774",
              dateStamp = "2018-08-25",
              name = "Harrison Specht, Guillaume Harmange, David H. Perlman, Edward Emmott, Zachary Niziolek, Bogdan Budnik, Nikolai Slavov",
              lab = "Slavov Lab",
              instrumentModel = "Q-Exactive Orbitrap",
              instrumentManufacturer = "Thermo Scientific",
              softwareName = "MaxQuant",
              softwareVersion = "1.6.2.3",
              switchingCriteria = "the top 5 most intense precursors with charges 2 to 4 were selected for HCD fragmentation at resolution 70,000 with a max fill time of 300ms.",
              ionSource = "not specified",
              ionSourceDetails = "not specified",
              analyser = "orbitrap",
              analyserDetails = "a 0.7 Th isolation window was used for MS2 scans.",
              collisionEnergy = "not specified")

# Create the MSnSet object
x <- new("MSnSet", exprs = edat, 
         featureData = fdat,
         phenoData = pdat,
         experimentData = expdat)

# Save data as Rda file
stopifnot(validObject(x))
assign("specht2018", x)
save(specht2018, file = file.path("../../data/specht2018.rda"),
     compress = "xz", compression_level = 9)
# Note: saving is assumed to occur in "(...)/scpdata/inst/scripts"



