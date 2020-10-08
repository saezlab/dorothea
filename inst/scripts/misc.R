library(biomaRt)
library(dplyr)
library(tibble)
library(readr)
library(tidyr)
library(readxl)
library(janitor)

# load tf census from Lambert et al (Document S1; Table S2) https://doi.org/10.1016/j.cell.2018.01.029
# download file
download.file("https://ars.els-cdn.com/content/image/1-s2.0-S0092867418301065-mmc2.xlsx", 
              destfile = "inst/extdata/annotations/lambert_tf_census.xlsx")
tf_census = read_excel("inst/extdata/annotations/lambert_tf_census.xlsx", 
                       sheet = 2, skip = 1) %>%
  clean_names() %>%
  rename(is_tf = x4, tf = name) %>%
  filter(is_tf == "Yes") %>%
  distinct(tf)
# remove downloaded file
file.remove("inst/extdata/annotations/lambert_tf_census.xlsx")

# load tf annotation from Garcia-Alonso et al (Table S1) http://www.genome.org/cgi/doi/10.1101/gr.240663.118
# file is downloaded manually from
tf_annotation = read_excel("inst/extdata/annotations/GarciaAlonso_supplemental_table_S1.xlsx",
                           skip = 1) %>%
  select(tf = TF, class = mode_of_regulation) %>%
  na_if("-")

# update tf annotation with census from lambert, mode of regulation/class is assumed to be NA for tfs not available yet in tf annotation
updated_tf_annotation = tf_census %>%
  left_join(tf_annotation, by="tf") 

saveRDS(updated_tf_annotation, "inst/extdata/annotations/tf_annotation.rds")


# Gene annotation with protein ids
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl", 
                  "http://aug2020.archive.ensembl.org")
biomart_output = getBM(mart = ensembl,
                       attributes=c("hgnc_symbol","uniprotswissprot",
                                    "ensembl_gene_id", "uniprot_gn_id"),
                       filters = "transcript_biotype",
                       values = c("protein_coding")) %>%

  as_tibble() %>%
  na_if("")

saveRDS(biomart_output, "inst/extdata/annotations/gene_annotation.rds")

# hgnc - mgi mapping
mouse_ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl", 
                        host = "http://aug2020.archive.ensembl.org")
human_ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl", 
                        host = "http://aug2020.archive.ensembl.org")
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
uniprot_anno = read_delim('inst/extdata/annotations/hsa_uniprot_20180314.txt',
                          delim = "\t") %>%
  janitor::clean_names() %>%
  drop_na(gene_names) %>%
  separate_rows(gene_names, sep="[ ;]") %>%
  group_by(entry, entry_name, ensembl_transcript) %>%
  mutate(gene_name = head(gene_names, 1)) %>%
  ungroup() %>%
  distinct(entry, entry_name, gene_name, ensembl_transcript)

saveRDS(uniprot_anno, "inst/extdata/annotations/uniprot_annotation.rds")

