library(dplyr)
library(purrr)
library(tidyr)
library(readr)
library(stringr)
library(tibble)

source("inst/scripts/utils.R")

tfs = load_tf_census()

anno = readRDS("inst/extdata/annotations/uniprot_annotation.rds") %>%
  distinct(entry, entry_name, gene_name)

hocomoco = read_delim("inst/extdata/tf_target_sources/tfbs/hocomoco/scanning_output/fimo/fimo_annotated.txt", delim = "\t") %>%
  janitor::clean_names() %>%
  select(-sequence_name, -start, -stop) %>%
  distinct()

jaspar = read_delim("inst/extdata/tf_target_sources/tfbs/jaspar/scanning_output/fimo/fimo_annotated.txt", delim = "\t") %>%
  janitor::clean_names() %>%
  select(-sequence_name, -start, -stop) %>%
  distinct()




# map hocomocos "number motif id" with "gene names" by uniprots "entry name"
mapping = hocomoco %>%
  distinct(number_motif_id) %>%
  # mutate(entry_name = str_split(number_motif_id, pattern = "\\.")[[1]][1]) %>%
  separate(number_motif_id, into = "entry_name", extra = "drop", remove = FALSE,
           sep = "\\.") %>%
  left_join(anno, by=c("entry_name")) %>%
  mutate(gene_name = case_when(is.na(gene_name) ~ str_remove(entry_name, "_HUMAN"),
                               TRUE ~ gene_name)) %>%
  dplyr::select(number_motif_id, entry_name, entry, gene_name) %>%
  distinct(number_motif_id, gene_name)

# integrate the mapping in the original hocomoco database
hocomoco_updated = hocomoco %>%
  select(-motif_alt_id) %>%
  left_join(mapping, by="number_motif_id") %>%
  rename(motif_alt_id = gene_name) %>%
  select(number_motif_id, motif_alt_id, everything())

databases = bind_rows(
  mutate(jaspar, source = "jaspar"),
  mutate(hocomoco_updated, source = "hocomoco")
  )



cutoff = 500
df = databases %>%
  group_by(source, number_motif_id) %>%
  slice(1:cutoff) %>%
  ungroup() %>%
  distinct(source, tf = motif_alt_id, target = gene) %>%
  mutate(tf = str_to_upper(tf)) %>%
  filter(tf %in% tfs) %>%
  drop_na() %>%
  arrange(tf, target)

df %>%
  filter(source == "jaspar") %>%
  distinct(tf, target) %>%
  saveRDS("inst/extdata/networks/tfbs/jaspar/network.rds")

df %>%
  filter(source == "hocomoco") %>%
  distinct(tf, target) %>%
  saveRDS("inst/extdata/networks/tfbs/hocomoco/network.rds")
