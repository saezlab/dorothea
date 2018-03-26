# DoRothEA (Discriminant Regulon Expression Analysis)


DoRothEA (Discriminant Regulon Expression Analysis) is a framework to estimate single sample TF activities from consensus TF-target DNA binding networks.  Please see [Garcia-Alonso et al 2018](http://cancerres.aacrjournals.org/content/early/2017/12/09/0008-5472.CAN-17-1679) for more datils.


## Usage

The main function is SLEA (Sample Level Enrichment Analisys), located in ``src/lib_enrichment_scores.r``: 


```
SLEA(E, genesets, method, M = NULL, permutations = 1000, filter_E = F)

```


```
  Info:
  - E = Expression matrix (rows=genes; columns=samples) scaled and recentered. Note that genes need to be in a comparable scale (e.g. z-transformed). 
  - genesets = list with two elemens: 1) NAME = vector with the names of the genes set; 2) GENES = named list of vectors containing the gene sets. Values indicate the edges (1 in this case) while names indicate the gene symbol of the targets. 
  - method = name of the method used for the SLEA ("MEAN"; "GSEAlm"; "VIPER"; "ssGSEA"; "GSVA").
  - M = Methylation binary matrix (optional). If null will be ignored
  - permutations. Number of permutations needed for the method "MEAN"
  - filter_E. Logical indicating if genes in E not in the gene set should be removed/ Default FALSE.
  - min_geneset_size. Minimum number of targets per TF. Default (n = 3).
  - max_shared. Maximum number of TFs per target. Default (n = 10).  
```

The SLEA function calls the following methods:
   - "MEAN" (aka z-SCORE). Calls a [z-test](https://genomemedicine.biomedcentral.com/articles/10.1186/gm327) method.
   - "GSEAlm" (aka MLR). Uses a [Multiple Linear Regression](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1182396/) approach.
   - "VIPER". Calls the [aREA](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5040167/) (analytic Rank-based Enrichment Analysis) method from the [VIPER](https://www.bioconductor.org/packages/release/bioc/html/viper.html) R package.
  - "ssGSEA". Calls the [ssGSEA](https://www.ncbi.nlm.nih.gov/pubmed/19847166)
     (single-sample GSEA) method from the [GSVA](https://bioconductor.org/packages/release/bioc/html/GSVA.html) R package. 
  - "GSVA". Calls the [GSVA](https://www.ncbi.nlm.nih.gov/pubmed/23323831)
       (Gene Set Variation Analysis) method from the [GSVA](https://bioconductor.org/packages/release/bioc/html/GSVA.html) R package.
       

See ``src/example.r`` for an example.

To reproduce the TF activities from [Garcia-Alonso et al 2018](https://www.ncbi.nlm.nih.gov/pubmed/29229604), please use the cl_voom_batchcor_dupmerged_KDCFgenenorm.rdata file in www.synapse.org [syn10463688](https://www.synapse.org/#!Synapse:syn10463688/wiki/463140). This file contains gene-wise normalized expression estimates for protein coding genes in more than 1,300 cell lines. The consensus TF regulons are available in the _data > regulons_ folder.


## About the TF regulons

DoRothEA makes use of consensus TF regulon (CTFR), composed of TF–target regulatory interactions from 13 public resources covering different TF-binding evidences, including TF-binding site (TFBS) predictions, chromatin immunoprecipitation coupled with high-throughput data (ChIP-X), text-mining derived and manually curated TF–target interactions. For each TF, the CTFR is defined by selecting TF–target interactions reported in more than one source. These TF–target interactions are unsigned and unweighted.

The regulons are provided as an R object at ``data/regulons/CTFRs_v122016.rdata``


## Citation

[Garcia-Alonso et al 2018](https://www.ncbi.nlm.nih.gov/pubmed/29229604)
Transcription Factor Activities Enhance Markers of Drug Sensitivity in Cancer.
Cancer Res February 1 2018 (78) (3) 769-780; 
DOI: 10.1158/0008-5472.CAN-17-1679


```
@article{garcia2018transcription,
  title={Transcription factor activities enhance markers of drug sensitivity in cancer},
  author={Garcia-Alonso, Luz and Iorio, Francesco and Matchan, Angela and Fonseca, Nuno and Jaaks, Patricia and Peat, Gareth and Pignatelli, Miguel and Falcone, Fiammetta and Benes, Cyril H and Dunham, Ian and others},
  journal={Cancer research},
  volume={78},
  number={3},
  pages={769--780},
  year={2018},
  publisher={American Association for Cancer Research}
}
```



