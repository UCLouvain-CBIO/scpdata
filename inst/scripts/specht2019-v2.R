
####---- Specht et al. 2019 ---####


# Specht, Harrison, Edward Emmott, Toni Koller, and Nikolai Slavov. 2019. 
# “High-Throughput Single-Cell Proteomics Quantifies the Emergence of Macrophage 
# Heterogeneity.” bioRxiv. https://doi.org/10.1101/665307.

# NOTE: the script is version 2 because the data set was updated by the authors

library(dplyr)
library(magrittr)
library(MSnbase)
library(tidyr)
setwd("./inst/scripts")

# The data was downloaded from https://scope2.slavovlab.net/docs/data to 
# scpdata/inst/extdata/specht2019-v2

# The end of the script contains some utility functions that must be run.


####---- Create the experiment data ----####

# Create the experiment data for the MSnSet object
# Information was found in the article
expdat <- new("MIAPE",
              title = "High-throughput single-cell proteomics quantifies the emergence of macrophage heterogeneity",
              abstract = "The fate and physiology of individual cells are controlled by protein interactions. Yet, our abilityto quantitatively analyze proteins in single cells has remained limited. To overcome this barrier, wedeveloped SCoPE2. It lowers cost and hands-on time by introducing automated and miniaturizedsample preparation while substantially increasing quantitative accuracy.  These advances enabledus to analyze the emergence of cellular heterogeneity as homogeneous monocytes differentiatedinto macrophage-like cells in the absence of polarizing cytokines.  SCoPE2 quantified over 2,700proteins in 1,018 single monocytes and macrophages in ten days of instrument time, and the quan-tified proteins allowed us to discern single cells by cell type.  Furthermore, the data uncovered acontinuous gradient of proteome states for the macrophage-like cells, suggesting that macrophageheterogeneity may emerge even in the absence of polarizing cytokines.  Parallel measurements oftranscripts by 10x Genomics scRNA-seq suggest that SCoPE2 samples 20-fold more copies pergene, thus supporting quantification with improved count statistics.  Joint analysis of the data indi-cated that most genes had similar responses at the protein and RNA levels, though the responses ofhundreds of genes differed.  Our methodology lays the foundation for automated and quantitativesingle-cell analysis of proteins by mass-spectrometry.",
              url = "https://doi.org/10.1101/665307",
              dateStamp = "2019-06-09",
              name = "Harrison Specht, Edward Emmott, Aleksandra A. Petelski, R. Gray Huffman, David H. Perlman, Marco Serra, Peter Kharchenko, Antonius Koller, Nikolai Slavov",
              lab = "Slavov Lab",
              instrumentModel = "Q-Exactive mass spectrometer",
              instrumentManufacturer = "Thermo Scientific",
              softwareName = "MaxQuant",
              softwareVersion = "1.6.2.3",
              switchingCriteria = "After a precursor scan from 450 to 1600 m/z at 70,000 resolving power, the top 7 most intense precursor ions with charges 2 to 4 and above the AGC min threshold of 20,000 were isolated for MS2 analysis via a 0.7 Th isolation window",
              ionSource = "ESI",
              ionSourceDetails = "Electrospray voltage was set to 2,200V, applied at the end of the analytical column. To reduce atmospheric background ions and enhance peptide signal to noise ratio, an Active Background Ion Reduction Device (ABIRD, by ESI Source Solutons, LLC, Woburn MA, USA) was used at the nanospray interface. The temperature of ion transfer tube was 250 degrees Celsius and the S-lens RF level set to 80.",
              analyser = "orbitrap",
              analyserDetails = "Precursor ions were accumulated for at most 300ms. Then they were fragmented via HCD at a and the fragments analyzed at 70,000 resolving power. Dynamic exclusion was used with a duration of 30 seconds with a mass tolerance of 10ppm.",
              collisionEnergy = "33 eV (normalized to m/z 500, z=1)")


####---- Load the data ----####

# MaxQuant output data
mqFile <- "../extdata/specht2019-v2/MSV000084660/quant/ev_SCoPE2_runs.csv"
sc0 <- readMSnSet2(file = mqFile, sep = ",", header = TRUE,
                   ecol = grepEcols(mqFile, split = ",",
                                    pattern = "intensity[.]\\d"))
sc <- sc0

# Load the sample annotations
annot <- read.csv("../extdata/specht2019-v2/MSV000084660/metadata/annotation.csv", 
                  row.names = 1, check.names = F)
batch <- read.csv("../extdata/specht2019-v2/MSV000084660/metadata/batch.csv", 
                  row.names = 1, check.names = F)

# Add a peptide sequence-charge field
fData(sc)$sequence_charge <- paste0(fData(sc)$Modified.sequence, 
                                    fData(sc)$Charge)


####---- Filter data based on sample annotation ----####

# Keep only sc runs and runs that
# * are not blanks, QC, reference, ladder, 10 cells, 100 cells, 1000 cells
# * are not in the run FP95, FP96, FP103
sel <- !grepl("blank|col19|col2[0-4]|arrier|QC|Ref|Master|SQC|blank|FP9[56]|FP103", 
              fData(sc)$Raw.file)
# * are not described in the metadata
sel <- sel & fData(sc)$Raw.file %in% rownames(annot)
# Remove sample not described in the metadata
sc.runs <- unique(fData(sc)$Raw.file[sel])
sc <- sc[sel, ]


####---- Filter peptides and runs ----####

# Remove the reverse hits (from decoy database) and contaminants
sc <- scp_filterPep(sc, thresPIF = 0.8)

# Remove runs with insufficient identified peptides 
sc <- scp_filterRun(sc, thresPep = 300)

# Remove peptides that are more than 10% the intensity of the carrier
sc <- scp_filterSCR(sc, samples = 4:16, carrier = 1, thresh = 0.1) 

# Remove spectra with poor identification confidence
# The PEP and q-values were updated using DART-ID
sc <- scp_filterFDR(sc, FDR = 0.01)


####---- Numeric normalization ----####

# Normalize single cell runs to normalization channel
sc <- scp_normalizeRI(sc, ref_col = 2)


####---- Format data to a peptide x sample matrix ----####

.keep <- c("Raw.file", "sequence_charge", "Modified.sequence", "Sequence", 
           "Length", "Proteins", "Leading.razor.protein", "Gene.names", 
           "Protein.names", "Mass")

# Convert data to long format 
sc <- scp_exprsToLong(sc, keepFdat = .keep)
# Remove duplicates
sc <- sc[!duplicated(fData(sc)[, c("sequence_charge", "Raw.file", "channel")]), ]
# Convert data to wide format
sc <- scp_exprsToWide(sc)
# dims = 15699 x 2816

####---- Missing data cleaning ----####

# Replace zero's by NA
sc <- scp_cleanMissing(sc, misVal = NA)


####---- Add sample metadata ----####

pd <- data.frame(row.names = sampleNames(sc))
pd$run <- sapply(rownames(pd), function(x) strsplit(x, "-")[[1]][1])
pd$channel <- sapply(rownames(pd), function(x) 
  gsub("Reporter[.]intensity[.]", "RI", strsplit(x, "-")[[1]][2])
)
pd$samp_type <- sapply(1:nrow(pd), function(i) annot[pd$run[i], pd$channel[i]])
pd <- cbind(pd, batch[pd$run, ])
pData(sc) <- pd


####---- Keep only single cell data ----####

sel <- pData(sc_nn)$samp_type %in% c("sc_u","sc_m0", "sc_0") &
  pData(sc_nn)$run %in% sc.runs
sc <- sc[sel, ]


####--- Filter cells ----####

sc <- scp_filterCells(sc_nn, sc_nc, sc_nc_lb, 
                       rowNorm = TRUE, npep = 6, 
                       qprobs = 0.3, q_thres = -2.5,
                       cv_thres = 0.43, 
                       median_thres = -1.3,
                       Plot = TRUE)
# Note1: I was not able to reproduce the CV calculation but had something very 
# similar
# Note2: npep or rowNorm has no big impact on the CV distribution. The 
# distribution changes a lot whether we use the carrier normalized or the 
# reference normalized data 
hist(exprs(sc_nnf), breaks = "FD", xlim = c(0, 2), col = "darkseagreen")


####---- Normalize rows and columns ----####

sc_nnfn <- scp_normalize_stat(sc_nnf, what = "column", fun = "/", stats = median)
sc_nnfn <- scp_normalize_stat(sc_nnfn, what = "row", fun = "/", stats = mean)

hist(exprs(sc_nnfn), breaks = "FD", xlim = c(-2, 2), col = "darkseagreen")


####---- Filter based on missingness ----####

sc_nnfn <- scp_filterNA(sc_nnfn, "row", pNA = 0.99)
sc_nnfn <- scp_filterNA(sc_nnfn, "column", pNA = 0.99)

hist(exprs(sc_nnfn), breaks = "FD", xlim = c(-2, 2), col = "darkseagreen")


####---- Log2 transform ----####

sc_final <- log(sc_nnfn, base = 2)

hist(exprs(sc_final), breaks = "FD", xlim = c(-2, 2), col = "darkseagreen")

# This sc_final is the peptide data as provided in Specth et al article.


####---- Utility functions ----####

scp_filterPep <- function(obj, thresPIF = 0.8){
  # Check arguments
  if(!inherits(obj, "MSnSet")) stop("'obj' must be an MSnSet object")
  # Filter out contamination peptides
  sc <- sc[is.na(fData(sc)$Reverse) | 
             fData(sc)$Reverse != "+", ]
  sc <- sc[!grepl("^REV|^CON", fData(sc)$Leading.razor.protein), ]
  sc <- sc[is.na(fData(sc)$Potential.contaminant) | 
             fData(sc)$Potential.contaminant != "+", ]
  # Filter out non pure spectra
  sc <- sc[!is.na(fData(sc)$PIF) & fData(sc)$PIF > thresPIF, ]
  # Log the process to the MSnSet object
  obj@processingData@processing <-
    c(processingData(obj)@processing,
      paste0("Remove contamination peptides and contaminated spectra (PIF <", 
             round(thresPIF, 3), "): ", date()))
  return(obj)
}

scp_filterRun <- function(obj, thresPep = 300){
  # Check arguments
  if(!inherits(obj, "MSnSet")) stop("'obj' must be an MSnSet object")
  # Filter out runs with less than thresPep peptides
  pep.t <- table(fData(sc)$Raw.file)
  sel <- fData(sc)$Raw.file %in% names(pep.t[pep.t >= 300])
  sc <- sc[sel, ]
  # Log the process to the MSnSet object
  obj@processingData@processing <-
    c(processingData(obj)@processing,
      paste0("Remove runs containing less than ", thresPep, "peptides (", 
             sum(!sel), " runs were removed): ", date()))
  return(obj)
}

# Filter peptdes based on the peptide false discovery rate (FDR)
scp_filterFDR <- function(obj, FDR = 0.01){
  # Check arguments
  if(!inherits(obj, "MSnSet")) stop("'obj' must be an MSnSet object")
  # Compute the FDR per sequence-charge 
  fdat <- fData(obj)
  calc_fdr <- function(pep)
    return((cumsum(pep[order(pep)]) / seq_along(pep))[order(order(pep))])
  fdat %<>% group_by(sequence_charge) %>% 
    mutate(pep_fdr = calc_fdr(dart_PEP))
  fData(obj) <- data.frame(fdat, row.names = featureNames(obj))
  # Filter peptides with high FDR
  obj <- obj[fData(obj)$pep_fdr < 0.01, ]
  # Log process to the MSnSet
  obj@processingData@processing <-
    c(processingData(obj)@processing,
      paste0("Filter peptides based on the FDR (threshold = ", round(FDR, 3), 
             "): ", date()))
  
  return(obj)
}

# Filter peptdes based on the sample over carrier ratio (SCR)
scp_filterSCR <- function(obj, samples, carrier, thresh = 0.1){
  # Check arguments
  if(!inherits(obj, "MSnSet")) stop("'obj' must be an MSnSet object")
  if(length(carrier) != 1) stop("'carrier' must have length 1. ")
  # Compute ratios
  ratios <- apply(exprs(obj)[, samples, drop = FALSE], 2, 
                  function(x) x/exprs(obj[, carrier]))
  
  # Filter data
  sel <- rowMeans(ratios, na.rm = TRUE) < thresh
  sel[is.na(sel)] <- FALSE
  obj <- obj[sel, ]
  # Log filtering to the MSnSet
  obj@processingData@processing <-
    c(processingData(obj)@processing,
      paste0("Remove features with a sample to carrier ratio higher than ",
             round(thresh * 100, 1), " %: ", date()))
  return(obj)
}


# Normalize to a reference channel
scp_normalizeRI <- scp_normaliseRI <- function(obj, ref_col){
  # Check arguments
  if(!inherits(obj, "MSnSet")) stop("'obj' must be an MSnSet object")
  if(length(ref_col) != 1) stop("'ref_col' must have length 1. ")
  if(is.numeric(ref_col)) ref_col <- sampleNames(obj)[ref_col]
  # Normalize channels
  exprs(obj) <- exprs(obj)/exprs(obj)[, ref_col]
  # Log the process to the MSnSet object
  obj@processingData@processing <-
    c(processingData(obj)@processing,
      paste0("Normalize reporter intensities by dividing every channel ",
             "by the '", ref_col, "' channel: ", date()))
  return(obj)
}


scp_exprsToLong <- function(obj, keepFdat = c("Raw.file", "sequence_charge")){
  # Check arguments
  if(!inherits(obj, "MSnSet")) stop("'obj' must be an MSnSet object")
  if(is.numeric((keepFdat))) keepFdat <- varLabels(featureData(obj))[keepFdat]
  if(!all(keepFdat %in% varLabels(featureData(obj)))) 
    stop(paste0("Wrong 'keepFdat' specification. The variable(s) '", 
                keepFdat[!keepFdat %in% varLabels(featureData(obj))], 
                "' is/are not present in 'varLabels(featureData(obj))'"))
  if(ncol(pData(obj)) != 0) stop("No implementation for elongating 'pData(obj)'. It should be a 0-row data frame.")
  # Combine expression and metadata
  ed <- exprs(obj) 
  fd <- fData(obj)[, keepFdat]
  dat <- cbind(fd, ed)
  # Elongate data
  dat <- pivot_longer(data = as.data.frame(dat), names_to = "channel", 
                      values_to = "Intensity", cols = colnames(ed))
  # Unbind expression and feature data
  ed <- dat$Intensity
  fd <- as.data.frame(dat[, colnames(dat) != "Intensity"])
  # Log process 
  obj@processingData@processing <- 
    c(processingData(obj)@processing, paste0("Convert data to long format: ", date()))
  # Create a new MSnSet object with updated information
  ObjNew <- new("MSnSet",
                featureData = AnnotatedDataFrame(fd),
                experimentData = experimentData(obj),
                processingData = processingData(obj),
                exprs = as.matrix(ed))
  stopifnot(validObject(ObjNew))
  return(ObjNew)
}

scp_exprsToWide <- function(obj){
  require(tidyr)
  # Check arguments
  if(!inherits(obj, "MSnSet")) stop("'obj' must be an MSnSet object")
  # Get data to widen
  ed <- exprs(obj) 
  fd <- fData(obj)
  # Create a unique id for every sample (will become the new column names)
  id <- paste0(fd$Raw.file, "-", fd$channel)
  # Widen expression data
  ed <- data.frame(seq = fd$sequence_charge, id = id, Intensity = as.vector(ed))
  ed <- pivot_wider(ed, id_cols = "seq", names_from = id,
                    values_from = "Intensity", 
                    values_fill = list("Intensity" = NA))
  rown <- ed[,1]
  ed <- as.matrix(ed[,-1])
  rownames(ed) <- rown$seq
  # Adapt feature data
  fd <- fd[!duplicated(fd$sequence_charge), ]
  rownames(fd) <- fd$sequence_charge
  fd <- fd[rownames(ed),]
  # Log process 
  obj@processingData@processing <- 
    c(processingData(obj)@processing, paste0("Convert data to wide format: ", date()))
  # Create a new MSnSet object with updated information
  ObjNew <- new("MSnSet",
                featureData = AnnotatedDataFrame(fd),
                experimentData = experimentData(obj),
                processingData = processingData(obj),
                exprs = ed)
  stopifnot(validObject(ObjNew))
  return(ObjNew)
}

scp_cleanMissing <- function(obj, misVal){
  # Check arguments
  if(!inherits(obj, "MSnSet")) stop("'obj' must be an MSnSet object")
  # Clean zero and infinite data points
  exprs(obj)[exprs(obj) == 0 | is.infinite(exprs(obj))] <- misVal
  # Log process 
  obj@processingData@processing <- 
    c(processingData(obj)@processing, paste0("Replace infinite values and 0's by ", 
                                             misVal, ": ", date()))
  return(obj)
}



scp_filterCells <- function(obj_nn,
                            obj_nc, 
                            obj_lb,
                            qprobs = 0.3,
                            q_thres = -2.5,
                            median_thres = -1.3,
                            cv_thres = 0.43,
                            rowNorm = TRUE, npep = 6, Plot = TRUE){
  # Check arguments
  # TODO 
  # Compute the quantile rRI per sample
  rRI_q <- apply(log10(exprs(obj_lb)), 2, quantile, 
                 probs = qprobs, na.rm = TRUE)
  rRI_q[is.infinite(rRI_q)] <- -3
  # Compute the median rRI per sample
  rRI_median <- apply(log10(exprs(obj_lb)), 2, median, na.rm = TRUE)
  rRI_median[is.infinite(rRI_median)] <- -3
  # Compute the median protein CV per sample
  cv_median <- computeCV(obj_nc, rowNorm = rowNorm, npep = npep)
  # Plot distributions
  if(Plot){
    par(mfrow = c(3,1))
    hist(rRI_q, breaks = 30, xlab = "Quantile rRI", col = "bisque3",
         main = paste0("Histogram of the quantile rRI (p = ", 
                       round(qprobs, 2), ")"))
    polygon(x = rep(c(q_thres, -100), each = 2), y = rep(c(-100, 100), 2),
            col = rgb(1, 0.2, 0.2, 0.2))
    hist(rRI_median, breaks = 30, main = "Histogram of the median rRI",
         xlab = "Median rRI", col = "bisque3")
    polygon(x = rep(c(median_thres, 100), each = 2), y = rep(c(-100, 100), 2),
            col = rgb(1, 0.2, 0.2, 0.2))
    hist(cv_median, breaks = 30, main = "Histogram of the median protein CV",
         xlab = "Median CV",  col = "bisque3")
    polygon(x = rep(c(cv_thres, 100), each = 2), y = rep(c(-100, 100), 2),
            col = rgb(1, 0.2, 0.2, 0.2))
  }
  # Filter data 
  sel <- cv_median < cv_thres & rRI_q > q_thres & rRI_median < median_thres
  sum(sel, na.rm = TRUE)
  obj_nn <- obj_nn[, sel]
  # Log the process
  obj_nn@processingData@processing <- 
    c(processingData(obj_nn)@processing, 
      paste0("Filter samples based on the log10 quantile rRI (p =", qprobs, 
             "; threshold = ", q_thres, "), the log10 median rRI (threshold = ", 
             median_thres, "), and the median protein CV (threshold = ", 
             cv_thres, "): ", date()))
  return(obj_nn)
}


