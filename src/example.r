

# Load functions
source('src/lib_enrichment_scores.r')

# Load TF regulon genesets
load('data/CTFRs_v122016.rdata')
# Load expression matrix
load('data/example_expression_mat.rdata')

# Estimate TF activities
TF_activities = SLEA(E = E, genesets = CTFRs_genesets, method = 'GSVA')$NES
