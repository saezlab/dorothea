library(dplyr)
library(purrr)
library(tidyr)
library(readr)
library(stringr)
library(tibble)

source("inst/scripts/utils.R")

tfs = load_tf_census()

#### GTEx ####
# Load networks
gtex_path = "inst/extdata/tf_target_sources/inferred/gtex/tissue_specific"


gtex_df = list.files(path = gtex_path, full.names = TRUE, recursive = TRUE,
                pattern = "viperRegulon") %>%
  map_dfr(function(path) {
    tissue = str_split(path, pattern = "/") %>%
      pluck(1,7)

    message(tissue)
    regulon = get(load(path))

    map_dfr(regulon, function(i) {
      tf_target = i$tfmode %>%
        enframe(name = "target", value = "mor") %>%
        mutate(tissue = tissue)
    }, .id = "tf") %>%
      mutate(tf = str_remove(tf, str_c(" ", tissue)))
  })

gtex_network = gtex_df %>%
  filter(tf %in% tfs) %>%
  mutate(mor = case_when(mor > 0 ~ 1,
                         mor < 0 ~ -1)) %>%
  # keep only interactions that are reported in at least in 3 tissues
  add_count(tf, target, mor, name = "signed_evidence") %>%
  filter(signed_evidence >= 3) %>%
  distinct(tf, target, mor) %>%
  # keep only interactions with unambigious mor
  add_count(tf, target) %>%
  filter(n == 1) %>%
  select(-n) %>%
  arrange(tf, target)

saveRDS(gtex_network,
        "inst/extdata/networks/inferred/gtex/pantissue/network.rds")

#### tcga ####
# Load networks
tcga_path = "inst/extdata/tf_target_sources/inferred/tcga/cancer_specific"


tcga_df = list.files(path = tcga_path, full.names = TRUE, recursive = TRUE,
                     pattern = "viperRegulon") %>%
  map_dfr(function(path) {
    cancer = basename(path) %>%
      str_split(pattern = "_") %>%
      pluck(1,1)

    message(cancer)
    regulon = get(load(path))

    map_dfr(regulon, function(i) {
      tf_target = i$tfmode %>%
        enframe(name = "target", value = "mor") %>%
        mutate(cancer = cancer)
    }, .id = "tf") %>%
      mutate(tf = str_remove(tf, str_c(" ","tcga_", cancer)))
  })

tcga_network = tcga_df %>%
  filter(tf %in% tfs) %>%
  mutate(mor = case_when(mor > 0 ~ 1,
                         mor < 0 ~ -1)) %>%
  # keep only interactions that are reported at least in 3 cancer types
  count(tf, target, mor, name = "signed_evidence") %>%
  filter(signed_evidence >= 3) %>%
  # keep only interactions with unambigious mor
  add_count(tf, target) %>%
  filter(n == 1) %>%
  dplyr::select(-n, -signed_evidence) %>%
  arrange(tf, target) %>%
  distinct(tf, target, mor)

saveRDS(tcga_network,
        "inst/extdata/networks/inferred/tcga/pancancer/network.rds")
