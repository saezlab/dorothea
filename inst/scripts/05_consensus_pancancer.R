library(dplyr)
library(purrr)
library(tidyr)
library(readr)
library(stringr)
library(tibble)
library(limma)

source("inst/scripts/utils.R")

tfs = load_tf_census()

#### Load networks from all types of evidences ####
n = list.files("inst/extdata/networks", pattern = "network",
               recursive = TRUE, full.names = TRUE) %>%
  # remove networks inferred from gtex data
  discard(.p = ~str_detect(.x, "gtex")) %>%
  map_dfr(function(path) {
    message(path)
    net = readRDS(path)

    # extract evidence type
    evidence = path %>%
      str_split("/") %>%
      pluck(1,4)

    # extract database
    database = path %>%
      str_split("/") %>%
      pluck(1,5)

    # set mode of regulation to 0 for interactions without any mor information
    if (!("mor" %in% colnames(net))) {
      net = net %>% mutate(mor = 0)
    }
    net = net %>%
      transmute(tf, target, mor,
                evidence = evidence,
                database = database)
    return(net)
  }) %>%
  distinct()

# update deprecated gene symbols of tfs and targets
target_aliases = n %>%
  distinct(target) %>%
  mutate(alias = alias2SymbolTable(target)) %>%
  mutate(alias = coalesce(alias, target))

tf_aliases = n %>%
  distinct(tf) %>%
  mutate(alias = alias2SymbolTable(tf)) %>%
  mutate(alias = coalesce(alias, tf))

n = n %>%
  inner_join(target_aliases, by="target") %>%
  rename(old_target = target) %>%
  select(tf, target = alias, mor, evidence, database, old_target) %>%
  inner_join(tf_aliases, by="tf") %>%
  rename(old_tf = tf) %>%
  select(tf = alias, target, mor, evidence, database, old_tf, old_target) %>%
  select(-old_tf, -old_target) %>%
  distinct()

#### Confidence class A ####
# interactions in >= 2 curated databases
a1 = n %>%
  filter(evidence == "curated") %>%
  count(tf, target) %>%
  filter(n >= 2) %>%
  select(tf, target)

# interactions in > 0 reviews/TF_e
a2 = n %>%
  filter(database %in% c("reviews", "tf_e")) %>%
  select(tf, target)

# signed interactions in curated and any other evidence (non curated)
signed = n %>%
  filter(evidence == "curated" & mor != 0) %>%
  distinct(tf, target)

noncurated = n %>%
  filter(evidence != "curated") %>%
  distinct(tf, target)

a3 = noncurated %>%
  inner_join(signed, by=c("tf", "target")) %>%
  distinct(tf, target)

# interactions in all 4 evidences
a4 = n %>%
  distinct(tf, target, evidence) %>%
  count(tf, target) %>%
  filter(n == 4) %>%
  select(-n)

a = bind_rows(a1, a2, a3, a4) %>%
  distinct() %>%
  mutate(confidence = "A")

##### Confidence class B ####
# interactions in curated databases and ChIP-seq
b1 = n %>%
  distinct(tf, target, evidence) %>%
  filter(evidence %in% c("curated", "chip_seq")) %>%
  count(tf, target) %>%
  filter(n==2) %>%
  select(-n)

# interaction in curated databases, predictions and TFBS
b2 = n %>%
  distinct(tf, target, evidence) %>%
  filter(evidence %in% c("curated", "tfbs", "inferred")) %>%
  count(tf, target) %>%
  filter(n==3) %>%
  select(-n)

# interaction in ChIP-seq, predictions and TFBS
b3 = n %>%
  distinct(tf, target, evidence) %>%
  filter(evidence %in% c("chip_seq", "tfbs", "inferred")) %>%
  count(tf, target) %>%
  filter(n==3) %>%
  select(-n)

b = bind_rows(b1, b2, b3) %>%
  distinct() %>%
  mutate(confidence = "B")


#### Confidence class C ####
# interaction in curated and TFBS
c1 = n %>%
  distinct(tf, target, evidence) %>%
  filter(evidence %in% c("curated", "tfbs")) %>%
  count(tf, target) %>%
  filter(n==2) %>%
  select(-n)

# interactions in ChIP_Seq and TFBS
c2 = n %>%
  distinct(tf, target, evidence) %>%
  filter(evidence %in% c("chip_seq", "tfbs")) %>%
  count(tf, target) %>%
  filter(n==2) %>%
  select(-n)

# interactions in ChIP_seq & inferred
c3 = n %>%
  distinct(tf, target, evidence) %>%
  filter(evidence %in% c("chip_seq", "inferred")) %>%
  count(tf, target) %>%
  filter(n==2) %>%
  select(-n)

c = bind_rows(c1, c2, c3) %>%
  distinct() %>%
  mutate(confidence = "C")

#### Confidence class D ####
# interactions occurring only in curated
d1 = n %>%
  distinct(tf, target, evidence) %>%
  filter(evidence == "curated") %>%
  distinct(tf, target)

# interactions occurring only in ChIP-seq
d2 = n %>%
  distinct(tf, target, evidence) %>%
  filter(evidence == "chip_seq") %>%
  distinct(tf, target)

d = bind_rows(d1, d2) %>%
  distinct() %>%
  mutate(confidence = "D")

#### Confidence class E ####
# interaction occuring only TFBS
e1 = n %>%
  distinct(tf, target, evidence) %>%
  filter(evidence == "tfbs") %>%
  distinct() %>%
  select(-evidence)

# interaction occuring only in inferred
e2 = n %>%
  distinct(tf, target, evidence) %>%
  filter(evidence == "inferred") %>%
  distinct() %>%
  select(-evidence)

e = bind_rows(e1, e2) %>%
  mutate(confidence = "E")

#### Combination of all confidence classes ####
# combine all interactions in decreasing confidence
# This leads to duplicated interactions, e.g. an interaction reported in two
# curated databases will have confidence level A (see a1) and D (see d1). We
# prioritize those interactions based on their confidence level (from A to E)
scored_interactions = bind_rows(a, b,c,d,e) %>%
  group_by(tf, target) %>%
  slice(1) %>%
  ungroup() %>%
  arrange(tf, target)


#### Fix sign conflict ####
# identify interactions in curated databases with contradictory mor (-1 vs 1)
interaction_with_sign_conflict = n %>%
  filter(evidence == "curated" & mor != 0) %>%
  distinct(tf, target, mor) %>%
  count(tf, target) %>%
  filter(n==2) %>%
  select(-n)

# do this interaction have further mor evidence in inferred/co-expression data
with_coexpression = n %>%
  filter(evidence == "inferred") %>%
  semi_join(interaction_with_sign_conflict, by=c("tf", "target"))


# if not check whether this contradictory interactions have a mor in curated
# databases with given priority among the databases
curated_database_priority = c("tf_e","tf_act", "trrust", "trrd_via_tf_act",
                              "nfi_regulome_db", "oreganno")

wo_coexpression = n %>%
  filter(evidence == "curated") %>%
  anti_join(with_coexpression, by=c("tf", "target")) %>%
  semi_join(interaction_with_sign_conflict, by=c("tf", "target")) %>%
  mutate(database = factor(database, levels = curated_database_priority)) %>%
  drop_na(database) %>%
  arrange(database, tf, target) %>%
  group_by(tf, target) %>%
  slice(1) %>%
  ungroup()

# interaction with changed mor
changed_sign = bind_rows(with_coexpression, wo_coexpression) %>%
  transmute(tf, target, mor, evidence = "curated")

# integrate the changed interaction into database
interactions = n %>%
  filter(tf %in% tfs) %>%
  left_join(changed_sign, by=c("tf", "target", "evidence")) %>%
  mutate(mor = case_when(mor.x == 0 ~ 0,
                            is.na(mor.y) ~ mor.x,
                            TRUE ~ mor.y)) %>%
  select(tf, target, mor, evidence, database) %>%
  mutate(mor = case_when(evidence == "inferred" ~ 0,
                            evidence != "inferred" ~ mor))

#### Build final database ####
# some interaction occur multiple times. Prioritize duplicated interactions
# based on mode of regulation and then on evidence type (see order of factor
# levels)
unique_interactions = interactions %>%
  mutate(evidence = factor(evidence, levels = c("curated", "chip_seq",
                                                "tfbs", "inferred")),
         mor = factor(mor, levels = c(1,-1, 0))) %>%
  arrange(mor, evidence) %>%
  group_by(tf, target) %>%
  slice(1) %>%
  ungroup() %>%
  select(tf, target, mor, evidence, database) %>%
  mutate(mor = as.numeric(as.character(mor)))

# retrieve database information for each interaction
which_databases = interactions %>%
  group_by(tf, target, evidence) %>%
  summarise(x = str_c(database, collapse = ",")) %>%
  ungroup() %>%
  mutate(key = case_when(
    evidence == "chip_seq" ~ "which_chip_seq",
    evidence == "curated" ~ "which_curated",
    evidence == "inferred" ~ "which_inferred",
    evidence == "tfbs" ~ "which_tfbs")) %>%
  select(-evidence) %>%
  spread(key, x, fill = "none")

# retrieve evidence information for each interaction
which_evidence = interactions %>%
  distinct(tf, target, evidence) %>%
  mutate(evidence = case_when(
    evidence == "chip_seq" ~ "is_evidence_chip_seq",
    evidence == "curated" ~ "is_evidence_curated",
    evidence == "inferred" ~ "is_evidence_inferred",
    evidence == "tfbs" ~ "is_evidence_tfbs"
  )) %>%
  mutate(val = TRUE) %>%
  spread(evidence, val, fill=FALSE)

# retrieve pubmed ids for interaction from curated databases
pubmed = list.files("inst/extdata/networks/curated", pattern = "pubmed",
                    full.names = TRUE, recursive = TRUE) %>%
  map_dfr(readRDS) %>%
  na_if("N,A")

pubmed_ids = pubmed %>%
  drop_na(pubmed_id) %>%
  separate_rows(pubmed_id, sep = ",") %>%
  distinct(tf, target, pubmed_id) %>%
  group_by(tf, target) %>%
  summarise(pubmed_id = str_c(pubmed_id, collapse = ",")) %>%
  ungroup()

 # combine all information (scores, which evidence,which database, pubmed)
entire_database_pancancer = unique_interactions %>%
  inner_join(scored_interactions, by = c("tf", "target")) %>%
  inner_join(which_evidence, by=c("tf", "target")) %>%
  inner_join(which_databases, by=c("tf", "target")) %>%
  left_join(pubmed_ids, by=c("tf", "target")) %>%
  replace_na(list(pubmed_id = "-")) %>%
  arrange(tf, target) %>%
  select(-evidence, -database) %>%
  # remove interactions exclusively reported by kegg
  filter(!(which_curated == "kegg" & which_chip_seq == "none" &
             which_inferred == "none" & which_tfbs == "none"))

# make final database with adapted mor. assign to all interactions mor of 1
# beside known repressors. Those get a mor of -1
tf_annotation = readRDS("inst/extdata/annotations/tf_annotation.rds") %>%
  filter(class == "repressors")
final_database_human_pancancer = entire_database_pancancer %>%
  distinct(tf, target, mor, confidence) %>%
  left_join(tf_annotation, by="tf") %>%
  mutate(mor = case_when(mor == 0 & is.na(class) ~ 1,
                         mor == 0 & class == "repressors" ~ -1,
                         TRUE ~ mor)) %>%
  select(-class)

#### Translate final database to mouse symbols ####
anno = readRDS("inst/extdata/annotations/hgnc_mgi_annotation.rds") %>%
  filter(!str_detect(mgi_symbol, "^Gm[:digit:]+"))

final_database_mouse_pancancer = final_database_human_pancancer %>%
  rename(hgnc_symbol = tf) %>%
  inner_join(anno, by="hgnc_symbol") %>%
  select(tf = mgi_symbol, target, mor, confidence) %>%
  # now translate targets
  rename(hgnc_symbol = target) %>%
  inner_join(anno, by="hgnc_symbol") %>%
  distinct(tf, mor, confidence, target = mgi_symbol) %>%
  group_by(tf, target) %>%
  # due to the mapping it can happen that now an interactions has several
  # confidence levels/mors. To be more conservative the lowest level/max
  # letter is chosen
  filter(confidence == max(confidence)) %>%
  # when multiple human gene symbols map to a single mouse symbol there can 
  # be ambiguous mor's (e.g. tf == "Egr1" & target == "Klk1b5"). For those cases
  # we assign a +1 as mor
  slice_max(order_by = mor, n=1) %>%
  ungroup() %>%
  distinct()

#### Top 10 database - Human ####
# To provide the most confident regulon for each TF, we aggregate the TF–target
# interactions with the highest possible confidence score that resulted in a
# regulon size equal to or greater than 10 targets.
top_10_database_human_pancancer = final_database_human_pancancer %>%
  add_count(tf, name = "total_interactions") %>%
  # filter out tfs with less than 4 targets
  filter(total_interactions >= 4) %>%
  # remove self-regulation of a TF
  filter(tf != target) %>%
  nest(regulon = -tf) %>%
  mutate(r = regulon %>% map(function(regulon) {
    c = regulon %>%
      arrange(confidence) %>%
      count(confidence) %>%
      mutate(cs = cumsum(n)) %>%
      filter(cs>=10)

    if (nrow(c) == 0) {
      summary_confidence = regulon %>%
        count(confidence) %>%
        slice_max(order_by = n, n = 1, with_ties = FALSE) %>%
        pull(confidence)

      res = regulon %>%
        transmute(summary_confidence, target, mor)

    } else {
      res = c %>%
        slice(1) %>%
        mutate(summary_confidence = confidence) %>%
        select(-n) %>%
        mutate(confidence = str_c(LETTERS[1:which(LETTERS == confidence)],
                                  collapse = ",")) %>%
        ungroup() %>%
        separate_rows(confidence, sep = ",") %>%
        inner_join(regulon, by="confidence") %>%
        select(summary_confidence, target, mor)
    }

    # if mor is unknown (0) assign a positive mode of regulation (1)
    res %>%
      mutate(mor = case_when(mor == 0 ~ 1,
                             TRUE ~ mor))
  })) %>%
  unnest(r) %>%
  select(-regulon) %>%
  rename(confidence = summary_confidence)

#### Top 10 database - Mouse ####
# To provide the most confident regulon for each TF, we aggregate the TF–target
# interactions with the highest possible confidence score that resulted in a
# regulon size equal to or greater than 10 targets.
top_10_database_mouse_pancancer = final_database_mouse_pancancer %>%
  add_count(tf, name = "total_interactions") %>%
  # filter out tfs with less than 4 targets
  filter(total_interactions >= 4) %>%
  # remove self-regulation of a TF
  filter(tf != target) %>%
  nest(regulon = -tf) %>%
  mutate(r = regulon %>% map(function(regulon) {
    c = regulon %>%
      arrange(confidence) %>%
      count(confidence) %>%
      mutate(cs = cumsum(n)) %>%
      filter(cs>=10)

    if (nrow(c) == 0) {
      summary_confidence = regulon %>%
        count(confidence) %>%
        slice_max(order_by = n, n = 1, with_ties = FALSE) %>%
        pull(confidence)

      res = regulon %>%
        transmute(summary_confidence, target, mor)

    } else {
      res = c %>%
        slice(1) %>%
        mutate(summary_confidence = confidence) %>%
        select(-n) %>%
        mutate(confidence = str_c(LETTERS[1:which(LETTERS == confidence)],
                                  collapse = ",")) %>%
        ungroup() %>%
        separate_rows(confidence, sep = ",") %>%
        inner_join(regulon, by="confidence") %>%
        select(summary_confidence, target, mor)
    }
    # if mor is unknown (0) assign a positive mode of regulation (1)
    res %>%
      mutate(mor = case_when(mor == 0 ~ 1,
                             TRUE ~ mor))
  })) %>%
  unnest(r) %>%
  select(-regulon) %>%
  rename(confidence = summary_confidence)

#### Save data ####
### Human
# top 10 database
dorothea_hs_pancancer = top_10_database_human_pancancer
save(dorothea_hs_pancancer, 
     file = "data/dorothea_hs_pancancer.rda", compress="xz")

### Mouse
# top 10 database
dorothea_mm_pancancer = top_10_database_mouse_pancancer
save(dorothea_mm_pancancer, 
     file = "data/dorothea_mm_pancancer.rda", compress="xz")
