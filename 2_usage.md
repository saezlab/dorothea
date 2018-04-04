---
layout: page
title: Usage
---

## DoRothEA (v1)
Below is an example how to run the first version of DoRothEA. The corresponding GitHub repository you find [here](https://github.com/saezlab/DoRothEA/releases/tag/version1)
```
# Load functions
source('src/lib_enrichment_scores.r')

# Load TF regulon genesets
load('data/CTFRs_v122016.rdata')
# Load expression matrix
load('data/example_expression_mat.rdata')

# Estimate TF activities
TF_activities = SLEA(E = E, genesets = CTFRs_genesets, method = 'GSVA')$NES
```
## DoRothEA (v2)
Below you find two examples how to calculate TF activities using the second version of DoRothEA.

```
# Load requeired packages
require(viper) 

# Load TF regulon genesets in VIPER format
load('data/TFregulons/Robjects_VIPERformat/consensus/BEST_viperRegulon.rdata')
# Clean TF names & explore object
names(viper_regulon) = sapply(strsplit(names(viper_regulon), split = ' - '), head, 1)
# Explore the regulons object
names(viper_regulon)[1:10]
viper_regulon[[1]]

##########################################################################################
## Example 1: Computing single-sample TF activities from a normalized gene expression matrix 
##########################################################################################
# Load expression matrix: NOTE that genes at comparable scales (e.g. zscores)
load('data/expression/example_expressionMatrix_zscores.rdata')
# Explore the matrix
E[1:5, 1:5]
# Estimate TF activities
TF_activities = viper(eset = E, regulon = viper_regulon, nes = T, method = 'none', minsize = 4, eset.filter = F)
# Save results
write.csv(TF_activities, file = 'TFactivities_example1.csv')

##########################################################################################
## Example 2: Computing TF activity changes from a differential gene expression signature
##########################################################################################
# Load differential expression signature
load('data/expression/example_differentialExpression_results.rdata')
# Explore the signature
DEsignature[1:5, ]
# Exclude probes with unknown or duplicated gene symbol
DEsignature = subset(DEsignature, Symbol != "" )
DEsignature = subset(DEsignature, ! duplicated(Symbol))
# Estimatez-score values for the GES. Cheeck VIPER manual for details
myStatistics = matrix(DEsignature$logFC, dimnames = list(DEsignature$Symbol, 'logFC') )
myPvalue = matrix(DEsignature$P.Value, dimnames = list(DEsignature$Symbol, 'P.Value') )
mySignature = (qnorm(myPvalue/2, lower.tail = FALSE) * sign(myStatistics))[, 1]
mySignature = mySignature[order(mySignature, decreasing = T)]
# Estimate TF activities
mrs = msviper(ges = mySignature, regulon = viper_regulon, minsize = 4, ges.filter = F)
TF_activities = data.frame(Regulon = names(mrs$es$nes),
                           Size = mrs$es$size[ names(mrs$es$nes) ], 
                           NES = mrs$es$nes, 
                           p.value = mrs$es$p.value, 
                           FDR = p.adjust(mrs$es$p.value, method = 'fdr'))
TF_activities = TF_activities[ order(TF_activities$p.value), ]
# Save results
write.csv(TF_activities, file = 'TFactivities_example2.csv')
