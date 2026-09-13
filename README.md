
# SpaMTP <img src="man/figures/logo.png" align="right" height="100" alt="" />

<!-- badges: start -->


## *R*-based User-Friendly Spatial Metabolomic, Transcriptomic, and Proteomic Data Analysis Tool

<br>

<!-- badges: end -->

SpaMTP is an R package for integrative analysis of spatial metabolomics and
spatial transcriptomics data. Raw and continuous MSI spectra remain in
Cardinal's `MSImagingArrays` (unaligned mass axes) or matter-backed
`MSImagingExperiment`; aligned or binned MSI uses
`SpatialExperiment`, including `rowData()` annotations, `colData()` pixel
metadata, `spatialCoords()`, `imgData()`, and paired transcriptomes in
`altExp()`. Seurat is an optional interoperability target rather than the
package's foundation.

SpaMTP provides (1) mass-to-charge ratio (m/z) metabolite annotation, (2)
downstream statistical analysis including differential metabolite expression
and pathway analysis, and (3) integrative spatial-omics analysis. Explicit
`asSpatialExperiment()`, `asCardinal()`, `seuratToSpatialExperiment()`, and
`spatialExperimentToSeurat()` entry points connect these workflows.

The analysis package does not require Seurat. Only explicit conversion
functions use the optional SeuratObject package; loading, preprocessing,
plotting and integration do not convert inputs automatically. See
[migration status](BIOCONDUCTOR_MIGRATION.md) for validation and scientific limits.

## Moving a Seurat workflow to Bioconductor

Convert once, then keep the result in a Bioconductor container:

```r
# xy has numeric x/y columns and row names matching Seurat pixel names.
spe <- seuratToSpatialExperiment(
    seuratObject, assay = "SPM", layer = "counts", coordinates = xy)
spe <- normalizeSMData(spe, "LogNormalize", verbose = FALSE)
spe <- scaleSMData(spe)
spe <- runMetabolicPCA(spe, slot = "logcounts")

# Non-spatial transcriptomics does not need invented coordinates.
sce <- seuratToSingleCellExperiment(seuratObject, assay = "RNA")
```

| Operation in a Seurat workflow | Native alternative |
|---|---|
| NormalizeData | normalizeSMData: normcounts/logcounts assays |
| ScaleData | scaleSMData: feature-wise centring/scaling in a scaled assay |
| RunPCA | runMetabolicPCA / scater::runPCA: reducedDim |
| VlnPlot | mzViolinPlot / scater::plotExpression |
| Cell/feature metadata | colData / rowData |
| Multiple modalities | Paired altExp entries |
| Subsetting cells | Standard brackets, or subsetSPM with colData/colLabels |
| FindMultiModalNeighbors | multiOmicIntegration supplies a **different**, equal-weight PCA embedding, not WNN |

Conversion preserves the requested exact layer, metadata, identities and
paired alternative assays. It does not transfer optical rasters, graphs or
reductions. Split Seurat v5 layers must be joined beforehand or selected
explicitly. Other assays without the same pixels are skipped with a warning.
The original layer name is retained: a converted data layer is still called
data, so either select it explicitly or normalize counts to create logcounts.
Scaling creates a dense matrix; it is not an out-of-memory operation.

Optical images can be attached with addSpatialImage(). To return to an
external Seurat workflow, explicitly call spatialExperimentToSeurat().
Historical Seurat/WNN tutorials belong to the published-workflow branch.

For runnable native examples, see [modern data access](vignettes/Modern_Data_Access.Rmd)
and the [full annotation pipeline](vignettes/Metabolite_Annotation_Pipeline.Rmd).
Spatial coordinates have one authoritative store: `spatialCoords()`.
Seurat conversion reads a single image/FOV through `GetTissueCoordinates()`
before considering historical metadata copies; multiple images require an
explicit selection. Native matrix and METASPACE imports no longer duplicate
coordinates in pixel metadata.

Please head to the [**SpaMTP website**](https://genomicsmachinelearning.github.io/SpaMTP/) for tutorials and documentation, including the [**RaMP 3.0 indexed metabolite annotation pipeline**](https://genomicsmachinelearning.github.io/SpaMTP/developmental/articles/Metabolite_Annotation_Pipeline.html). Active development and website deployment remain in the [**project development repository**](https://github.com/GenomicsMachineLearning/SpaMTP/tree/developmental); this repository contains the streamlined Bioconductor package source.

SpaMTP is now published in *Nature Methods*: [**SpaMTP: integrative statistical analysis and visualization of spatial metabolomics and transcriptomics data**](https://doi.org/10.1038/s41592-026-03140-8).

<br>

<img src="man/figures/SpaMTP_Fig.png" alt="" style="background-color: white;" />

<br>

## Installation

Keep the coordinated source versions together: SpaMTP >= 0.99.5,
SpaMTPdb >= 0.99.4 and SpaMTPData >= 0.99.5. For local sibling checkouts,
install the data packages before the software package:

```sh
R CMD INSTALL ../SpaMTPdb
R CMD INSTALL ../SpaMTPData
R CMD INSTALL .
```

The [three-package workflow](vignettes/SpaMTP_Companion_Workflow.Rmd) runs
without downloads using SpaMTPData's native synthetic SpatialExperiment,
SpaMTPdb's registry and SpaMTP's analysis/annotation functions. Data packages
remain independent of the software package; conversion logic stays in SpaMTP.
The RaMP 3.0.7 snapshot and historical experiment-data 1.0.0 remain unchanged.
SpaMTPData now defaults to resource release 1.1.0, with seven native
SpatialExperiment RDS in [Zenodo record 22733262](https://zenodo.org/records/22733262)
and eleven unchanged auxiliary/Cardinal resources. Use `spaMTPData()` to
retrieve native datasets without Seurat. The one-time preparation recipe is
retained for reproducibility; historical archives require an explicit
`version = "1.0.0"` request.

For human gene mapping, see the [gene identifier workflow](vignettes/Gene_Identifier_Mapping.Rmd).
SpaMTPdb supplies a separately versioned, checksum-verified HGNC archive;
SpaMTPData supplies experiments with species provenance. SpaMTP's
`buildGeneMappingIndex()`, `mapGeneIdentifiers()` and
`annotateGeneIdentifiers()` resolve symbols and stable IDs with an input audit.
Pathway analyses unite memberships of RaMP records belonging to one HGNC gene
and count that gene once. Ambiguous aliases and conflicting gene records remain
explicitly unresolved. Targeted-panel Fisher analysis requires the full measured
panel as `universe`; all eligible pathways enter the BH correction.

During Bioconductor review, install the companion annotation package and this
submission source from GitHub:

``` r
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")

remotes::install_github("BCRL-tylu/SpaMTPdb")
remotes::install_github("BCRL-tylu/SpaMTPData")
remotes::install_github("SpaMTP-project/SpaMTP")
```

`SpaMTPdb` supplies versioned RaMP annotation and pathway resources.
`SpaMTPData`, installed separately when demonstration data are required,
provides ExperimentHub access to the larger objects used in tutorials.

After Bioconductor acceptance, the supported installation command will be:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("SpaMTP")
```

For tutorials and more information please visit the [SpaMTP website](https://genomicsmachinelearning.github.io/SpaMTP/)

## Citation

If SpaMTP contributes to your work, please cite:

> Causer, A., Lu, T., Kriel, J. *et al.* SpaMTP: integrative statistical analysis and visualization of spatial metabolomics and transcriptomics data. *Nature Methods* **23**, 1501–1506 (2026). https://doi.org/10.1038/s41592-026-03140-8

The citation can also be retrieved directly in R:

``` r
citation("SpaMTP")
```

## Maintainer

SpaMTP is currently maintained by **Tianyao Lu** ([GitHub](https://github.com/BCRL-tylu), [email](mailto:lu.t@wehi.edu.au)). **Andrew Causer** remains credited as an original author and former maintainer.
