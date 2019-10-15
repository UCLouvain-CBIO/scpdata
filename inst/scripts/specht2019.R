############################
#### Specht et al. 2019 ####
############################

# Specht, Harrison, Edward Emmott, Toni Koller, and Nikolai Slavov. 2019. 
# “High-Throughput Single-Cell Proteomics Quantifies the Emergence of Macrophage 
# Heterogeneity.” bioRxiv. https://doi.org/10.1101/665307.

library(MSnbase)

# Download data
dat <- read.csv(file = "http://slavovlab.net/scope2/data/Peptides-raw.csv")
protn <- dat$protein
pepn <- dat$peptide

# Create the MSnSet object holding the SCP data from Specht et al. 2019
x <- readMSnSet2(dat, 
                 ecol = -c(1,2)) # 1 and 2 are protein and peptide names, resp.
# TODO add cell data !!

# Add experiment data to the MSnSet object
experimentData(x) <- new("MIAPE",
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
                         analyser = "ion trap",
                         analyserDetails = "Precursor ions were accumulated for at most 300ms. Then they were fragmented via HCD at a and the fragments analyzed at 70,000 resolving power. Dynamic exclusion was used with a duration of 30 seconds with a mass tolerance of 10ppm.",
                         collisionEnergy = "33 eV (normalized to m/z 500, z=1)")

# Annotate fields
featureNames(x) <- fData(x)[, 1] 

# Save data as Rda file
stopifnot(validObject(x))
assign("specht2019", x)
save(specht2019, file = file.path("../../data/specht2019.rda"),
     compress = "xz", compression_level = 9)


