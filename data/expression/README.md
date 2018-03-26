## About the data

The purpose of the data in this folder is to provide gene expression signature examples to test the code.


   - example_expressionMatrix_zscores.rdata is an R object containing a z-tranformed gene expression matrix.

   - example_differentialExpression_results.rdata is an R object containing a differential expression signature derived using the R package limma. 
   
   
To reproduce the TF activities from [Garcia-Alonso et al 2018](https://www.ncbi.nlm.nih.gov/pubmed/29229604), please use the cl_voom_batchcor_dupmerged_KDCFgenenorm.rdata file in www.synapse.org [syn10463688](https://www.synapse.org/#!Synapse:syn10463688/wiki/463140). This file contains gene-wise normalized expression estimates for protein coding genes in more than 1,300 cell lines. The consensus TF regulons are available in the _data > regulons_ folder.