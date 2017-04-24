---
layout: page
title: Installation
---



SLAPenrich implements a statistical framework to identify pathways that tend to be recurrently genomically altered across a the samples of a genomic dataset. Differently from traditional over-recurrence analyses, SLAPenrich does not require the genes belonging to a given pathway to be statistically enriched among those altered in the individual samples. Consistently with the mutual exclusivity principle, and differently from other proposed computational tools, our approach assumes that the functionality of a given pathway might be altered in an individual sample if at least one of its genes is genomically altered. The method accounts for the differences in the **mutation rates** between samples and the **exonic lengths** of the genes in the pathways. It statistically tests against the null hypothesis that no associations between a pathway and the disease population under study does exist, assessing analytically the divergence of the total number of samples with alterations in a given pathway from its expectation. Moreover, the used formalism allows SLAPenrich to perform differential enrichment analysis of pathway alterations across different clinically relevant sub-populations of samples. SLAPenrich also includes function to visualise the identified enriched pathway implementing a heuristic sorting to highlight mutual exclusivity trends among the pattern of alterations of the composing genes.



## Installation

To install the `SLAPenrich` package from github, first run the following into your R console to install and load `devtools` package:

```r
install.packages("devtools")
library("devtools")
```

Next, install the and load `SLAPenrich` package from github:

```r
install_github("saezlab/SLAPenrich")
library(SLAPenrich)
```
