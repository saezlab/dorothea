library(dplyr)
library(purrr)
library(tidyr)
library(readr)
library(stringr)
library(tibble)

source("inst/scripts/utils.R")

tfs = load_tf_census()

remap_network = read_delim(
  "inst/extdata/tf_target_sources/chip_seq/remap/gene_tf_pairs_genesymbol.txt",
  delim = "\t", col_names = c("tf", "target", "enseml_target", "score")
  ) %>%
  filter(tf %in% tfs & score > 100)



c = 500
network = remap_network %>%
  group_by(tf) %>%
  top_n(c, score) %>%
  mutate(min_score = min(score)) %>%
  filter(score > min_score) %>%
  ungroup() %>%
  select(tf, target) %>%
  arrange(tf, target)

saveRDS(network, "inst/extdata/networks/chip_seq/remap/network.rds")
