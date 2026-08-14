
# SpaMTP <img src="man/figures/logo.png" align="right" height="100" alt="" />

<!-- badges: start -->


## *R*-based User-Friendly Spatial Metabolomic, Transcriptomic, and Proteomic Data Analysis Tool

<br>

<!-- badges: end -->

SpaMTP is an R package designed for the integrative analysis of spatial metabolomics and spatial transcriptomics data. SpaMTP inherits functionalities from two well established R packages (Cardinal and Seurat) to present a user-friendly platform for integrative spatial-omics analysis. Build on the foundation of a [*Seurat Class Object*](https://satijalab.org/seurat/), this package has three major functionalities which include; (1) mass-to-charge ratio (m/z) metabolite annotation, (2) various downstream statistical analysis including differential metabolite expression and pathway analysis, and (3) integrative spatial-omics analysis. In addition, this package includes various functions for data visualisation and data import/export, permitting flexible usage with other established R and Python  packages.   

Please head to the [**SpaMTP website**](https://genomicsmachinelearning.github.io/SpaMTP/) for tutorials and documentation, including the [**RaMP 3.0 indexed metabolite annotation pipeline**](https://genomicsmachinelearning.github.io/SpaMTP/developmental/articles/Metabolite_Annotation_Pipeline.html). Active development and website deployment remain in the [**project development repository**](https://github.com/GenomicsMachineLearning/SpaMTP/tree/developmental); this repository contains the streamlined Bioconductor package source.

SpaMTP is now published in *Nature Methods*: [**SpaMTP: integrative statistical analysis and visualization of spatial metabolomics and transcriptomics data**](https://doi.org/10.1038/s41592-026-03140-8).

<br>

<img src="man/figures/SpaMTP_Fig.png" alt="" style="background-color: white;" />

<br>

## Installation

During Bioconductor review, install the companion annotation package and this
submission source from GitHub:

``` r
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")

remotes::install_github("BCRL-tylu/SpaMTPdb")
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
