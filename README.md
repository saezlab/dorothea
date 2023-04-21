# DoRothEA: collection of human and mouse regulons <img src="man/figures/tool_logo.png" align="right" width="120" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/saezlab/dorothea/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/saezlab/dorothea/actions)
[![BioC status](http://bioconductor.org/shields/build/release/data-experiment/dorothea.svg)](https://bioconductor.org/checkResults/release/data-experiment-LATEST/dorothea)
[![codecov](https://codecov.io/gh/saezlab/dorothea/branch/master/graph/badge.svg)](https://codecov.io/gh/saezlab/dorothea)
![GitHub](https://img.shields.io/github/license/saezlab/dorothea)
<!-- badges: end -->

## Overview
DoRothEA is a gene regulatory network (GRN) containing signed transcription factor
(TF) - target gene interactions. DoRothEA regulons, the collection of a TF and
its transcriptional targets, were curated and collected from different types of
evidence for both human and mouse. A confidence level was assigned to each 
TF-target interaction based on the number of supporting evidence. 

<img src="man/figures/overview.png" align="center" width="600">

These regulons, coupled with any statistical method, can be
used to infer TF activities from bulk or single-cell transcriptomics. 

This is an R package for storing the regulons. To infer TF
activities, please check out
[decoupleR](https://doi.org/10.1093/bioadv/vbac016), available in
[R](https://saezlab.github.io/decoupleR/) or
[python](https://github.com/saezlab/decoupler-py).

## Installation

DoRothEA is available in
[Bioconductor](http://bioconductor.org/packages/release/data/experiment/html/dorothea.html). 
In addition, one can install the development version from the Github repository:
```r
## To install the package from Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("dorothea")

## To install the development version from the Github repo:
devtools::install_github("saezlab/dorothea")
```

## Updates

Since the original release, we have implemented some extensions in DoRothEA:

1. **Extension to mouse**:
  Originally DoRothEA was developed for the application to human data. 
  In a benchmark study we showed that DoRothEA is also applicable to mouse data, 
  as described in 
  [Holland et al., 2019](https://doi.org/10.1016/j.bbagrm.2019.194431). 
  Accordingly, we included new parameters to run mouse version of DoRothEA by 
  transforming the human genes to their mouse orthologs.
2. **Extension to single-cell RNA-seq data**:
  We showed that DoRothEA can be applied to scRNA-seq data, as described in
  [Holland et al., 2020](https://doi.org/10.1186/s13059-020-1949-z)

3. **Extension to other databases**
  We have released a new literature based GRN with increased coverage and better
  performance at identifying perturbed TFs, called [CollecTRI](https://github.com/saezlab/CollecTRI).
  We encourage users to use CollecTRI instead of DoRothEA. Vignettes on how to
  obtain activities are available at the [decoupleR package](https://saezlab.github.io/decoupleR/).

## License
DoRothEA is intended only for academic use as in contains resources
whose licenses don't permit commercial use. However, we developed a non-academic
version of DoRothEA by removing those resources (mainly TRED from the curated
databases). You can find the non-academic package with the regulons [here](https://github.com/saezlab/dorothea/tree/non-academic).

## Citation
If you use the DoRothEA resource, please cite:

> Garcia-Alonso L, Holland CH, Ibrahim MM, Turei D, Saez-Rodriguez J.
Benchmark and integration of resources for the estimation of human
transcription factor activities. _Genome Research._ 2019. DOI: [10.1101/gr.240663.118](https://doi.org/10.1101/gr.240663.118).

If you infer TF activities, please cite:

> Badia-i-Mompel P., Vélez Santiago J., Braunger J., Geiss C., Dimitrov
  D., Müller-Dott S., Taus P., Dugourd A., Holland C.H., Ramirez Flores R.O.
  and Saez-Rodriguez J. decoupleR: Ensemble of computational methods to
  infer biological activities from omics data. 2022. _Bioinformatics Advances_.
  DOI: [10.1093/bioadv/vbac016](https://doi.org/10.1093/bioadv/vbac016)

If you use the CollecTRI resource, please cite:

> Müller-Dott S., Tsirvouli E., Vázquez M., Ramirez Flores R.O.,
  Badia-i-Mompel P., Fallegger R., Lægreid A. and Saez-Rodriguez J.
  Expanding the coverage of regulons from high-confidence prior knowledge for
  accurate estimation of transcription factor activities. _bioRxiv_. 2023. DOI:
  [10.1101/2023.03.30.534849](https://doi.org/10.1101/2023.03.30.534849)
