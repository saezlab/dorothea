# dorothea <img src="man/figures/tool_logo.png" align="right" height="139">

<!-- badges: start -->
[![Build Status](https://travis-ci.com/saezlab/dorothea.svg?branch=master)](https://travis-ci.com/saezlab/dorothea)
[![codecov](https://codecov.io/gh/saezlab/dorothea/branch/master/graph/badge.svg)](https://codecov.io/gh/saezlab/dorothea)
![GitHub](https://img.shields.io/github/license/saezlab/dorothea)
<!-- badges: end -->

## Overview
DoRothEA is a gene set resource containing signed transcription factor (TF) - target interactions first described in [Garcia-Alonso et al., 2019](https://doi.org/10.1101/gr.240663.118). The collection of a TF and its transcriptional targets is defined as regulon. DoRothEA regulons were curated and collected from different types of evidence such as literature curated resources, ChIP-seq peaks, TF binding site motifs and interactions inferred directly from gene expression. 

<img src="man/figures/overview.png" align="center" width="400">

For each TF-target interaction we assigned a confidence level based on the number of supporting evidence. The confidence assigment comprises five levels, ranging from A (highest confidence) to E (lowest confidence). Interactions that are supported by all four lines of evidence, manually curated by experts in specific reviews, or supported both in at least two curated resources are considered to be highly reliable and were assigned an A level. Level B-D are reserved for curated and/or ChIP-seq interactions with different levels of additional evidence. Finally, E level is used for interactions that are uniquely supported by computational predictions. To provide the most confident regulon for each TF, we aggregated the TF-target interactions with the highest possible confidence score that resulted in a regulon size equal to or greater than ten targets. The final confidence level assigned to the TF regulon is the lowest confidence score of its component targets.

DoRothEA regulons can be coupled with several statistical method yielding a *functional analysis*-tool to infer TF activity from gene expression data. The activity is computed by considering not the gene expression of the TFs itself but the mRNA levels of their direct transcriptional targets. We define the transcriptional targets to as footprints of a TF on gene expression. A more detailed description of the concept of *footprint-based* analysis is available in the review [Dugourd et al., 2019](https://doi.org/10.1016/j.coisb.2019.04.002).

Typcially, DoRothEA is coupled with the statistical method [VIPER](https://www.bioconductor.org/packages/release/bioc/html/viper.html) as it incorporates the mode of regulation of each TF-target interaction. However, VIPER can be replaced by any other statistical method that aims to analyse gene sets, e.g. GSEA.


## Installation

```r
# install from bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("dorothea")

# install the development version from GitHub
# install.packages("devtools")
devtools::install_github("saezlab/dorothea")
```

## Updates
### Update #1 - extension to mouse
Originally DoRothEA contained only human regulons. In a benchmark study we showed that DoRothEA in combination with VIPER is also applicable to mouse data, as described in [Holland et al., 2019](https://doi.org/10.1016/j.bbagrm.2019.194431). Accordingly, we developed a mouse version of DoRothEA by transforming the human genes to their mouse orthologs.

### Update #2 - extension to single-cell RNA-seq data
Recent technological advances in single-cell RNA-seq enable the profiling of gene expression at the individual cell level. We showed that DoRothEA in combination with VIPER can be applied to scRNA-seq data, as described in [Holland et al., 2020](https://doi.org/10.1186/s13059-020-1949-z).

## Citing DoRothEA
Beside the original paper there are two additional papers expanding the usage of DoRothEA regulons.

* If you use DoRothEA for your research please cite the original publication: 
> Garcia-Alonso L, Holland CH, Ibrahim MM, Turei D, Saez-Rodriguez J. "Benchmark and integration of resources for the estimation of human transcription factor activities." _Genome Research._ 2019. DOI: [10.1101/gr.240663.118](https://doi.org/10.1101/gr.240663.118).

* If you use the mouse version of DoRothEA please cite additionally:
> Holland CH, Szalai B, Saez-Rodriguez J. "Transfer of regulatory knowledge from human to mouse for functional genomics analysis." _Biochimica et Biophysica Acta (BBA) - Gene Regulatory Mechanisms._ 2019. DOI: [10.1016/j.bbagrm.2019.194431](https://doi.org/10.1016/j.bbagrm.2019.194431).

* If you apply DoRothEA's regulons on single-cell RNA-seq data please cite additionally:
> Holland CH, Tanevski J, Perales-Patón J, Gleixner J, Kumar MP, Mereu E, Joughin BA, Stegle O, Lauffenburger DA, Heyn H, Szalai B, Saez-Rodriguez, J. "Robustness and applicability of transcription factor and pathway analysis tools on single-cell RNA-seq data." _Genome Biology._ 2020. DOI: [10.1186/s13059-020-1949-z](https://doi.org/10.1186/s13059-020-1949-z).


## Usage of DoRothEA for commercial purpose
DoRothEA as it stands is intended only for academic use as in contains resources whose licenses don't permit commerical use. However, we developed a non-academic version of DoRothEA by removing the critically resources (namely KEGG and TRED from the curated databases). You find the non-academic package with the regulons [here](https://github.com/saezlab/dorothea/tree/non-academic).
