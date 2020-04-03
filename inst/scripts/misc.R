library(biomaRt)
library(dplyr)
library(tibble)
library(readr)
library(tidyr)


# Gene annotation with protein ids
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
biomart_output = getBM(mart = ensembl,
                       attributes=c("hgnc_symbol","uniprotswissprot",
                                    "ensembl_gene_id", "uniprot_gn_id"),
                       filters = "transcript_biotype",
                       values = c("protein_coding")) %>%

  as_tibble() %>%
  na_if("")

saveRDS(biomart_output, "inst/extdata/annotations/gene_annotation.rds")

# hgnc - mgi mapping
mouse_ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")
human_ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# mgi symbol - hgnc symbol
hgnc_mgi_mapping = getLDS(attributes = c("mgi_symbol"),
       mart = mouse_ensembl,
       attributesL = c("hgnc_symbol"), martL = human_ensembl) %>%
  as_tibble() %>%
  rename(mgi_symbol = MGI.symbol, hgnc_symbol = HGNC.symbol) %>%
  na_if("") %>%
  drop_na()

saveRDS(hgnc_mgi_mapping, "inst/extdata/annotations/hgnc_mgi_annotation.rds")



# Uniprot entry annotation
uniprot_anno = read_delim('data/annotations/hsa_uniprot_20180314.txt',
                          delim = "\t") %>%
  janitor::clean_names() %>%
  drop_na(gene_names) %>%
  separate_rows(gene_names, sep="[ ;]") %>%
  group_by(entry, entry_name, ensembl_transcript) %>%
  mutate(gene_name = head(gene_names, 1)) %>%
  ungroup() %>%
  distinct(entry, entry_name, gene_name, ensembl_transcript)

saveRDS(uniprot_anno, "inst/extdata/annotations/uniprot_annotation.rds")

