# Single Cell Proteomics data sets

## The `scpdata` package

`scpdata` disseminates mass spectrometry (MS)-based single-cell
proteomics (SCP) data sets formatted using the `scp` data structure. The
data structure is `QFeatures` from the
[QFeatures](https://rformassspectrometry.github.io/QFeatures/), and is
described in the [`scp`
vignette](http://www.bioconductor.org/packages/release/bioc/vignettes/scp/inst/doc/scp.md).
The vignette describing [supported input formats for
`readQFeatures()`](https://rformassspectrometry.github.io/QFeatures/articles/readQFeatures2.html)
(equally applicable to single-cell proteomics data) describes how to
import quantitative data tables from farious different third-party
software.

In this vignette, we describe how to access the SCP data sets. To start,
we load the `scpdata` package.

[`library`](https://rdrr.io/r/base/library.html)`(`[`"scpdata"`](https://UCLouvain-CBIO.github.io/scpdata/)`)`

## Load data from `ExperimentHub`

The data is stored and managed through the
[`ExperimentHub`](https://bioconductor.org/packages/release/bioc/html/ExperimentHub.html)
infrastructure. We first create a connection with `ExperimentHub`.

`eh`` ``<-`` `[`ExperimentHub`](https://rdrr.io/pkg/ExperimentHub/man/ExperimentHub-class.html)`(``)`

You can list all data sets available in `scpdata` using the query
function.

`query``(``eh``, ``"scpdata"``)`` ``#> ExperimentHub with 32 records`` ``#> # snapshotDate(): 2026-09-04`` ``#> # $dataprovider: MassIVE, PRIDE, SlavovLab website, Dataverse, Zenodo`` ``#> # $species: Homo sapiens, Mus musculus, Rattus norvegicus, Gallus gallus`` ``#> # $rdataclass: QFeatures, SingleCellExperiment`` ``#> # additional mcols(): taxonomyid, genome, description,`` ``#> # coordinate_1_based, maintainer, rdatadateadded, preparerclass, tags,`` ``#> # rdatapath, sourceurl, sourcetype `` ``#> # retrieve records with, e.g., 'object[["EH3899"]]' `` ``#> `` ``#> title `` ``#> EH3899 | specht2019v2 `` ``#> EH3900 | specht2019v3 `` ``#> EH3901 | dou2019_lysates `` ``#> EH3902 | dou2019_mouse `` ``#> EH3903 | dou2019_boosting`` ``#> ... ... `` ``#> EH9610 | hu2023_K562 `` ``#> EH9611 | hu2023_oocyte `` ``#> EH9627 | ai2025a `` ``#> EH10516 | leduc2025 `` ``#> EH10517 | wu2026`

Another way to get information about the available data sets is to call
[`scpdata()`](https://UCLouvain-CBIO.github.io/scpdata/reference/scpdata.md).
This will retrieve all the available metadata. For example, we can
retrieve the data set titles along with the description to make an
informed choice about which data set to choose.

`info`` ``<-`` `[`scpdata`](https://UCLouvain-CBIO.github.io/scpdata/reference/scpdata.md)`(``)`` ``knitr``::`[`kable`](https://rdrr.io/pkg/knitr/man/kable.html)`(``info``[``, `[`c`](https://rdrr.io/r/base/c.html)`(``"title"``, ``"description"``)``]``)`

|  | title | description |
|:---|:---|:---|
| EH3899 | specht2019v2 | SCP expression data for monocytes (U-937) and macrophages at PSM, peptide and protein level |
| EH3900 | specht2019v3 | SCP expression data for more monocytes (U-937) and macrophages at PSM, peptide and protein level |
| EH3901 | dou2019_lysates | SCP expression data for Hela digests (0.2 or 10 ng) at PSM and protein level |
| EH3902 | dou2019_mouse | SCP expression data for C10, SVEC or Raw cells at PSM and protein level |
| EH3903 | dou2019_boosting | SCP expression data for C10, SVEC or Raw cells and 3 boosters (0, 5 or 50 ng) at PSM and protein level |
| EH3904 | zhu2018MCP | Near SCP expression data for micro-dissection rat brain samples (50, 100, or 200 µm width) at PSM level |
| EH3905 | zhu2018NC_hela | Near SCP expression data for HeLa samples (aproximately 12, 40, or 140 cells) at PSM level |
| EH3906 | zhu2018NC_lysates | Near SCP expression data for HeLa lysates (10, 40 and 140 cell equivalent) at PSM level |
| EH3907 | zhu2018NC_islets | Near SCP expression data for micro-dissected human pancreas samples (control patients or type 1 diabetes) at PSM level |
| EH3908 | cong2020AC | SCP expression data for Hela cells at PSM, peptide and protein level |
| EH3909 | zhu2019EL | SCP expression data for chicken utricle samples (1, 3, 5 or 20 cells) at PSM, peptide and protein level |
| EH6011 | liang2020_hela | Expression data for HeLa cells (0, 1, 10, 150, 500 cells) at PSM, peptide and protein level |
| EH7085 | schoof2021 | Single-cell proteomics data from OCI-AML8227 cell culture to reconstruct the cellular hierarchy. |
| EH7295 | williams2020_lfq | Single-cell label free proteomics data from a MCF10A cell line culture. |
| EH7296 | williams2020_tmt | Single-cell proteomics data from three acute myeloid leukemia cell line culture (MOLM-14, K562, CMK). |
| EH7712 | derks2022 | Single-cell and bulk (100-cell) proteomics data of PDAC, melanoma cells and monocytes. |
| EH7713 | brunner2022 | Single-cell proteomics data of cell cycle stages in HeLa. |
| EH8301 | leduc2022_pSCoPE | Single-cell proteomics data of 878 melanoma cells and 877 monocytes (pSCoPE). |
| EH8302 | leduc2022_plexDIA | Single-cell proteomics data of 126 melanoma cells (plexDIA). |
| EH8303 | woo2022_macrophage | Single-cell proteomics data from LPS-treated macrophages. |
| EH8304 | woo2022_lung | Single-cell proteomics data from primary human lung cells. |
| EH9450 | gregoire2023_mixCTRL | Single-cell proteomics data from two monocyte cell lines |
| EH9477 | khan2023 | Single-cell proteomics data of 421 MCF-10A cells undergoing EMT triggered by TGF-beta |
| EH9487 | guise2024 | Single-cell proteomics data of 108 postmortem CTL or ALS spinal moto neurons |
| EH9497 | petrosius2023_mES | Mouse embryonic stem cells across ground-state (m2i) and differentiation-permissive (m15) culture conditions. |
| EH9498 | petrosius2023_AstralAML | Single-cell proteomics data of 4 cell types from the OCI-AML8227 model. |
| EH9609 | krull2024 | Single-cell proteomics data IFN-Y response of U-2 OS cells |
| EH9610 | hu2023_K562 | Single-cell proteomics data of K562 cells |
| EH9611 | hu2023_oocyte | Single-cell proteomics data of oocytes |
| EH9627 | ai2025a | Single-cell proteomics data of adult cardiomyocytes from Ai et al. (2025) |
| EH10516 | leduc2025 | Principles of protein abundance regulation across single cells in a mammalian tissue Leduc et al, (2025) |
| EH10517 | wu2026 | Label-free single-cell proteomics data from prenatal human brain tissue (gestational weeks 13, 15, and 19) from Wu et al. (2026) |

To get one of the data sets (*e.g.* `dou2019_lysates`) you can either
retrieve it using the `ExperimentHub` query function

`scp`` ``<-`` ``eh``[[``"EH3901"``]``]`` ``#> see ?scpdata and browseVignettes('scpdata') for documentation`` ``#> loading from cache`` ``scp`` ``#> An instance of class QFeatures (type: scp) with 4 sets:`` ``#> `` ``#> [1] Hela_run_1: SingleCellExperiment with 24562 rows and 10 columns `` ``#> [2] Hela_run_2: SingleCellExperiment with 24310 rows and 10 columns `` ``#> [3] peptides: SingleCellExperiment with 13934 rows and 20 columns `` ``#> [4] proteins: SingleCellExperiment with 1641 rows and 20 columns`

or you can the use the built-in functions from `scpdata`

`scp`` ``<-`` `[`dou2019_lysates`](https://UCLouvain-CBIO.github.io/scpdata/reference/dou2019_lysates.md)`(``)`` ``#> see ?scpdata and browseVignettes('scpdata') for documentation`` ``#> loading from cache`` ``scp`` ``#> An instance of class QFeatures (type: scp) with 4 sets:`` ``#> `` ``#> [1] Hela_run_1: SingleCellExperiment with 24562 rows and 10 columns `` ``#> [2] Hela_run_2: SingleCellExperiment with 24310 rows and 10 columns `` ``#> [3] peptides: SingleCellExperiment with 13934 rows and 20 columns `` ``#> [4] proteins: SingleCellExperiment with 1641 rows and 20 columns`

## Data sets information

Each data set has been extensively documented in a separate man page
(*e.g.*
[`?dou2019_lysates`](https://UCLouvain-CBIO.github.io/scpdata/reference/dou2019_lysates.md)).
You can find information about the data content, the acquisition
protocol, the data collection procedure as well as the data sources and
reference.

## Data manipulation

For more information about manipulating the data sets, check the
[`scp`](http://bioconductor.org/packages/release/bioc/html/scp.md)
package. The `scp`
[vignette](http://bioconductor.org/packages/release/bioc/vignettes/scp/inst/doc/scp.md)
will guide you through a typical SCP data processing workflow. Once your
data is loaded from `scpdata` you can skip section 2 *Read in SCP data*
of the `scp` vignette.

## Session information

    R version 4.6.1 (2026-06-24)
    Platform: x86_64-pc-linux-gnu
    Running under: Ubuntu 24.04.4 LTS

    Matrix products: default
    BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0

    locale:
     [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
     [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
     [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
     [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
     [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

    time zone: UTC
    tzcode source: system (glibc)

    attached base packages:
    [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    [8] base     

    other attached packages:
     [1] scpdata_1.21.3              ExperimentHub_3.3.2        
     [3] AnnotationHub_4.3.2         BiocFileCache_3.3.0        
     [5] dbplyr_2.6.0                QFeatures_1.23.1           
     [7] MultiAssayExperiment_1.39.1 SummarizedExperiment_1.43.0
     [9] Biobase_2.73.2              GenomicRanges_1.65.4       
    [11] Seqinfo_1.3.2               IRanges_2.47.5             
    [13] S4Vectors_0.51.9            BiocGenerics_0.59.12       
    [15] generics_0.1.4              MatrixGenerics_1.25.0      
    [17] matrixStats_1.5.0           BiocStyle_2.41.0           

    loaded via a namespace (and not attached):
     [1] tidyselect_1.2.1            dplyr_1.2.1                
     [3] blob_1.3.0                  Biostrings_2.81.8          
     [5] filelock_1.0.3              SingleCellExperiment_1.35.2
     [7] fastmap_1.2.0               lazyeval_0.2.3             
     [9] digest_0.6.39               lifecycle_1.0.5            
    [11] cluster_2.1.8.3             ProtGenerics_1.45.0        
    [13] KEGGREST_1.53.6             RSQLite_3.53.3             
    [15] magrittr_2.0.5              compiler_4.6.1             
    [17] rlang_1.3.0                 sass_0.4.10                
    [19] tools_4.6.1                 igraph_2.3.3               
    [21] yaml_2.3.12                 knitr_1.51                 
    [23] S4Arrays_1.13.0             htmlwidgets_1.6.4          
    [25] bit_4.6.0                   curl_8.0.0                 
    [27] DelayedArray_0.39.6         plyr_1.8.9                 
    [29] abind_1.4-8                 withr_3.0.3                
    [31] purrr_1.2.2                 desc_1.4.3                 
    [33] grid_4.6.1                  MASS_7.3-66                
    [35] cli_3.6.6                   crayon_1.5.3               
    [37] rmarkdown_2.32              ragg_1.5.2                 
    [39] otel_0.2.0                  httr_1.4.9                 
    [41] reshape2_1.4.5              BiocBaseUtils_1.15.1       
    [43] DBI_1.3.0                   cachem_1.1.0               
    [45] stringr_1.6.0               AnnotationDbi_1.75.2       
    [47] AnnotationFilter_1.37.0     BiocManager_1.30.27        
    [49] XVector_0.53.0              vctrs_0.7.3                
    [51] Matrix_1.7-6                jsonlite_2.0.0             
    [53] bookdown_0.48               bit64_4.8.6                
    [55] clue_0.3-68                 systemfonts_1.3.2          
    [57] tidyr_1.3.2                 jquerylib_0.1.4            
    [59] glue_1.8.1                  pkgdown_2.2.1.9000         
    [61] stringi_1.8.9               BiocVersion_3.24.0         
    [63] tibble_3.3.1                pillar_1.11.1              
    [65] rappdirs_0.3.4              htmltools_0.5.9            
    [67] R6_2.6.1                    httr2_1.3.0                
    [69] textshaping_1.0.5           evaluate_1.0.5             
    [71] lattice_0.23-1              png_0.1-9                  
    [73] memoise_2.0.1               bslib_0.12.0               
    [75] Rcpp_1.1.2                  SparseArray_1.13.2         
    [77] xfun_0.60                   MsCoreUtils_1.25.4         
    [79] fs_2.1.0                    pkgconfig_2.0.3            

## License

This vignette is distributed under a [CC BY-SA
license](https://creativecommons.org/licenses/by-sa/2.0/).
