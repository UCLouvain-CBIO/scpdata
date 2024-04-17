##' Williams et al. 2020 (Anal. Chem.): 3 AML cell line
##'
##' Single-cell label data from three acute myeloid
##' leukemia cell line culture (MOLM-14, K562, CMK). The data were
##' acquired using a TMT-based quantification protocole and the
##' nanoPOTS technology. The objective was to demonstrate successful
##' use of the new nanoPOTS autosampler presented in the source
##' article. The samples contain either carrier (10 ng), reference
##' (0.2ng), empty or single-cell samples..
##'
##' @format A [QFeatures] object with 4 assays, each assay being a
##' [SingleCellExperiment] object:
##'
##' - `peptides_[intensity or corrected]`: 2 assays containing peptide
##'   reporter ion intensities or corrected reporter ion intensities
##'   as computed by MaxQuant.
##' - `proteins_[intensity or corrected]`: 2 assays containing protein
##'   reporter ion intensities or corrected reporter ion intensities
##'   as computed by MaxQuant.
##'
##' Sample annotation is stored in `colData(williams2020_tmt())`.
##'
##' @section Acquisition protocol:
##'
##' The data were acquired using the following setup. More information
##' can be found in the source article (see `References`).
##'
##' - **Sample isolation**: cultured MOLM-14, K562 or CMK cells were
##'   isolated using flow-cytometry based cell sorting and deposit on
##'   nanoPOTS microwells
##' - **Sample preparation**: cells are lysed using using a DDM lysis
##'   buffer. Proteins are digested with trypsin followed by TMT
##'   labelling and quanching with HA. The samples are then acidified
##'   with FA, pooled in a single samples (adding carrier and reference
##'   peptide mixtures), and dried until LC-MS/MS analysis.
##' - **Liquid chromatography**: peptides are loaded using the new
##'   autosampler described in the paper. Samples are loaded using a
##'   a homemade miniature syringe pump. The samples are then desalted
##'   and concentrated through a SPE column (4cm x 100µm i.d. packed
##'   with 5µm C18) with microflow LC pump. The peptides are then eluted
##'   from a long LC column (60cm x 50 µm i.d. packed with 3µm C18)
##'   coupled to a nanoflox LC pump at 150nL/mL (elution time is not
##'   expliceted).
##' - **Mass spectrometry**: MS/MS was performed on an Orbitrap Fusion
##'   Lumos Tribrid MS coupled to a 2kV ESI. MS1 setup: Orbitrap
##'   analyzer at resolution 120.000, AGC target of 1E6, accumulation
##'   of 246ms. MS2 setup: Orbitrap with HCD at resolution 120.000, AGC
##'   target of 1E6, accumulation of 246ms.
##' - **Raw data processing**: preprocessing using Maxquant v1.6.2.10
##'   that use Andromeda search engine (with UniProtKB 2016-21-29).
##'
##' @section Data collection:
##'
##' All data were collected from the MASSIVE repository (accession ID:
##' MSV000085230).
##'
##' The peptide and protein data were extracted from the
##' `Peptides_AML_SingleCell.txt` or `ProteinGroups_AML_SingleCell.txt`
##' files, respectively, in the `AML_SingleCell` folders.
##'
##' The tables were duplicated so that intensisities and corrected
##' intensities are contained in separate tables. Tables are then
##' converted to [SingleCellExperiment] objects. Sample annotations
##' were inferred from the sample names, from table S2 and from the
##' Experimental Section of the paper. All data is combined in
##' a [QFeatures] object. [AssayLinks] were stored between peptide
##' assays and their corresponding proteins assays based on the
##' leading razor protein (hence only unique peptides are linked to
##' proteins).
##'
##' The script to reproduce the `QFeatures` object is available at
##' `system.file("scripts", "make-data_williams2020_tmt.R", package = "scpdata")`
##'
##' @source
##'
##' The PSM and protein data can be downloaded from the MASSIVE
##' repository MSV000085230.
##'
##' @references
##'
##' **Source article**: Williams, Sarah M., Andrey V. Liyu, Chia-Feng
##' Tsai, Ronald J. Moore, Daniel J. Orton, William B. Chrisler,
##' Matthew J. Gaffrey, et al. 2020. “Automated Coupling of
##' Nanodroplet Sample Preparation with Liquid Chromatography-Mass
##' Spectrometry for High-Throughput Single-Cell Proteomics.”
##' Analytical Chemistry 92 (15): 10588–96.
##' ([link to article](http://dx.doi.org/10.1021/acs.analchem.0c01551)).
##'
##' @examples
##' \donttest{
##' williams2020_tmt()
##' }
##'
##' @keywords datasets
##'
"williams2020_tmt"
