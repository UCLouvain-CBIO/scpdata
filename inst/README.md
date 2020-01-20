
# Methodology used for collecting the data


1. Identify the data source and the annotations from the article, ask authors for data if needed. Add a summary of the available information in the section below
2. Create a new R script that converts the data to `SingleCellExperiment` objects
3. Add data documentation and the data collection procedure in `scpdata/R/data.R`
4. Update the `scpdata/README.md` file

# Available data sets in literature

### Dou et al 2019, Analytical Chemistry

- Software used: MS-GF+, [MASIC](https://github.com/PNNL-Comp-Mass-Spec/MASIC/releases/tag/v3.0.7111), R ( [RomicsProcessor](https://github.com/PNNL-Comp-Mass-Spec/RomicsProcessor), sva, pmartR, missMDA, FactoMineR ) 
- [Raw data](ftp://massive.ucsd.edu/MSV000084110/raw/)
- [Peptide data](ftp://massive.ucsd.edu/MSV000084110/other/)
- [Protein data](https://pubs.acs.org/doi/10.1021/acs.analchem.9b03349)
- Script: not available but used RomicsProcessor
- [Article](http://dx.doi.org/10.1021/acs.analchem.9b03349): Dou, Maowei, Geremy Clair, Chia-Feng Tsai, Kerui Xu, William B. Chrisler, Ryan L. Sontag, Rui Zhao, et al. 2019. “High-Throughput Single Cell Proteomics Enabled by Multiplex Isobaric Labeling in a Nanodroplet Sample Preparation Platform.” Analytical Chemistry 91 (20): 13119–27.
- In scpdata: **missing peptide data**

### Schoof et al 2019, BioRxiv

- Software used: Proteome, Discoverer, SCeptre, SCANPY
- [Raw data](https://www.ebi.ac.uk/pride/archive/projects/PXD015112/private): Username: reviewer87620@ebi.ac.uk; Password: 6BnVxQ9F. 
- [Peptide data](https://www.ebi.ac.uk/pride/archive/projects/PXD015112/private): Username: reviewer87620@ebi.ac.uk; Password: 6BnVxQ9F. 
- Protein data: from docker image kuikuisven/sceptre
- Scripts: [kuikuisven](https://github.com/kuikuisven)/sceptre (not available yet)
- [Article](http://dx.doi.org/10.1101/745679): Schoof, Erwin M., Nicolas Rapin, Simonas Savickas, Coline Gentil, Eric Lechman, James Seymour Haile, Ulrich auf Dem Keller, John E. Dick, and Bo T. Porse. 2019. “A Quantitative Single-Cell Proteomics Approach to Characterize an Acute Myeloid Leukemia Hierarchy.” bioRxiv. https://doi.org/10.1101/745679.
- In scpdata: **no**

### Zhu et al 2019, BioRxiv

- Software used: MaxQuant, Andromeda, R (limma, lmerTest, sva, CellTrails)
- [Raw data](ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/11/PXD014256)
- [Peptide data](ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/11/PXD014256): contained in the `SEARCH.zip`file
- [Protein data](ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/11/PXD014256): contained in the `OTHER.zip`file
- **(Scripts)**
- [Article](http://dx.doi.org/10.1101/727412): Zhu, Ying, Mirko Scheibinger, Daniel C. Ellwanger, Jocelyn F. Krey, Dongseok Choi, Ryan T. Kelly, Stefan Heller, and Peter G. Barr-Gillespie. 2019. “Single-Cell Proteomics Reveals Downregulation of TMSB4X to Drive Actin Release for Stereocilia Assembly.” bioRxiv. https://doi.org/10.1101/727412.
- In scpdata: **no**

### Specht et al 2019, BioRxiv

**Version 2**

- Software used: Maxquant, R (+DART-ID)
- [Raw data](ftp://massive.ucsd.edu/MSV000084660/raw/)
- [Peptide data](ftp://massive.ucsd.edu/MSV000084660/quant/)
- [Protein data](https://drive.google.com/open?id=1c5Z3b_2gOwDyHCLm9ycY3hXckY1GDd5L)
- [Scripts](https://github.com/SlavovLab/SCoPE2/tree/master/code). The script contains the full analysis from peptide to protein and generates all the figures from the article
- [Article](https://www.biorxiv.org/content/10.1101/647545v3.article-info): Emmott, Edward, Harrison Specht, Aleksandra Petelski, R. Gray Huffman, and Nikolai Slavov. 2019. “SCoPE2 for High-Throughput Single-Cell Quantitative Proteomics.” Preprint.
- In scpdata: **no**

**Version 1**

- Software used: Maxquant, R (+DART-ID)
- [Raw data](ftp://massive.ucsd.edu/MSV000083945/raw/scope2_raw/)
- [Peptide data](ftp://massive.ucsd.edu/MSV000083945/quant/), see also [here](https://drive.google.com/drive/folders/1cMQ-SIGpHwSfx9wJF2fIa-t8yX329LPM)
- [Protein data](https://scope2.slavovlab.net/docs/data)
- Scripts: version 1 is overwritten by version 2. The version 1 script can be obtained by calling the last September git [commit](https://github.com/SlavovLab/SCoPE2/commit/0f95bf4cb92df9a01168b25927f464d7194a5752) or from my [forked copy](https://github.com/UCLouvain-CBIO/scpScripts/tree/master/20191115-Specht2019/replicate%20specht/code). The script contains the full analysis from peptide to protein and generates all the figures from the article
- [Article](https://www.biorxiv.org/content/10.1101/665307v2): Emmott, Edward, Harrison Specht, Aleksandra Petelski, R. Gray Huffman, and Nikolai Slavov. 2019. “SCoPE2 for High-Throughput Single-Cell Quantitative Proteomics.” Preprint.
- In scpdata: yes

**scRNA-Seq**

scRNA-Seq [data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE142392) is also available (same data set in version 1 and 2)

### Specht et al 2018, BioRxiv

- Software used: Maxquant, R (+DART-ID)
- [Raw data](https://drive.google.com/drive/folders/15E8bkArJ5tV5gRJ3_o4DiowFJxzscVO_)
- [Peptide data](https://drive.google.com/drive/folders/19YG70I52DH5yETcZagdUjNZWNPs0JXVr?usp=sharing), but annotations are missing
- **(Protein data)**
- **(Scripts)**
- [Article](http://dx.doi.org/10.1101/399774): Specht, Harrison, Guillaume Harmange, David H. Perlman, Edward Emmott, Zachary Niziolek, Bogdan Budnik, and Nikolai Slavov. 2018. “Automated Sample Preparation for High-Throughput Single-Cell Proteomics.” bioRxiv. https://doi.org/10.1101/399774.
- In scpdata: yes

### Dou et al 2018, Chemical Science

- Sofwtare used: MaxQuant, Perseus
- [Raw data](ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2018/10/PXD010150)
- **(Peptide data)**
- **(Protein data)**
- **(Scripts)**
- [Article](https://pubs.rsc.org/en/content/articlehtml/2018/sc/c8sc02680g): Dou, Maowei, Ying Zhu, Andrey Liyu, Yiran Liang, Jing Chen, Paul D. Piehowski, Kerui Xu, et al. 2018. “Nanowell-Mediated Two-Dimensional Liquid Chromatography Enables Deep Proteome Profiling Of< 1000 Mammalian Cells.” Chemical Science  9 (34): 6944–51.
- In scpdata: **no**

### Budnik et al 2018, Genome Biology

- Software used: Maxquant, R
- [Raw data](ftp://massive.ucsd.edu/MSV000082077/raw/)
- [Peptide data](ftp://massive.ucsd.edu/MSV000082077/search/Search_Results/), but annotations are missing
- **(Protein data)**
- **(scripts)**
- [Article](http://dx.doi.org/10.1186/s13059-018-1547-5): Budnik, Bogdan, Ezra Levy, Guillaume Harmange, and Nikolai Slavov. 2018. “SCoPE-MS: Mass Spectrometry of Single Mammalian Cells Quantifies Proteome Heterogeneity during Cell Differentiation.” Genome Biology 19 (1): 161.
- In scpdata: **no**

### Zhu et al 2018, Angew Chem Int Ed Engl

- Software used: Maxquant, R 
- Raw data: [HeLa](ftp://proteomics@ftp.pnl.gov/outgoing/Laurent/Angew_Single_Cell/HeLa) and [Lung](ftp://proteomics@ftp.pnl.gov/outgoing/Laurent/Angew_Single_Cell/Lung_Cell). Username: proteomics; Password: Amt23Data
- **(Peptide data)**
- **(Protein data)**
- **(Scripts)**
- [Article](http://dx.doi.org/10.1002/anie.201802843): Zhu, Ying, Geremy Clair, William B. Chrisler, Yufeng Shen, Rui Zhao, Anil K. Shukla, Ronald J. Moore, et al. 2018. “Proteomic Analysis of Single Mammalian Cells Enabled by Microfluidic Nanodroplet Sample Preparation and Ultrasensitive NanoLC-MS.” Angewandte Chemie  57 (38): 12370–74.
- In scpdata: **no**

### Zhu et al 2018, Analytical Chemistry

- Software used: MaxQuant, Excel, OriginLab 2017
- [Raw data](ftp://proteomics@ftp.pnl.gov/outgoing/Laurent/OHSU_CTC). Username: proteomics; Password:  Amt23Data
- **(Peptide data)**
- **(Protein data)**
- **(Scripts)**
- [Article](http://dx.doi.org/10.1021/acs.analchem.8b03268): Zhu, Ying, Jennifer Podolak, Rui Zhao, Anil K. Shukla, Ronald J. Moore, George V. Thomas, and Ryan T. Kelly. 2018. “Proteome Profiling of 1 to 5 Spiked Circulating Tumor Cells Isolated from Whole Blood Using Immunodensity Enrichment, Laser Capture Microdissection, Nanodroplet Sample Processing, and Ultrasensitive nanoLC-MS.” Analytical Chemistry 90 (20): 11756–59.
- In scpdata: **no**

### Zhu et al 2018, Nature communication 

- Software used: MaxQuant, Perseus, OriginLab 2017
- [Raw data](ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2018/01/PXD006847): Cultured cells, islets, and Vail data are available
- [Peptide data](ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2018/01/PXD006847): Cultured cells, islets, and Vail data are available
- **(Protein data)**
- **(Scripts)**
- [Article](http://dx.doi.org/10.1038/s41467-018-03367-w): Zhu, Ying, Paul D. Piehowski, Rui Zhao, Jing Chen, Yufeng Shen, Ronald J. Moore, Anil K. Shukla, et al. 2018. “Nanodroplet Processing Platform for Deep and Quantitative Proteome Profiling of 10-100 Mammalian Cells.” Nature Communications 9 (1): 882.
- In scpdata: **yes**

### Zhu et al 2018, Mol Cell Proteomics

- Software used: MaxQuant, Andromeda, Perseus, OriginPro  2017
- [Raw data](ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2018/07/PXD008844)
- [Peptide data](ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2018/07/PXD008844)
- **Protein data**
- **(Scripts)**
- [Article](http://dx.doi.org/10.1074/mcp.TIR118.000686): Zhu, Ying, Maowei Dou, Paul D. Piehowski, Yiran Liang, Fangjun Wang, Rosalie K. Chu, William B. Chrisler, et al. 2018. “Spatially Resolved Proteome Mapping of Laser Capture Microdissected Tissue with Automated Sample Transfer to Nanodroplets.” Molecular & Cellular Proteomics: MCP 17 (9): 1864–74.
- In scpdata: **yes**

