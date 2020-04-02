library(KEGGgraph)
library(KEGGREST)
library(org.Hs.eg.db)
library(dplyr)
library(purrr)
library(tidyr)
library(readr)
library(stringr)
library(tibble)

source("inst/scripts/utils.R")

tfs = load_tf_census()

# litreature curated resources
#### Fantom4 ####
# load fantom table
fantom = read_delim(
  "inst/extdata/tf_target_sources/curated/fantom_4/edge.GoldStd_TF.tbl.txt",
  delim ="\t", skip = 24) %>%
  rename(pubmed_id = sym_value)

# load annotation table for fantom db
fantom_anno = read_delim(
  "inst/extdata/tf_target_sources/curated/fantom_4/feature.entrez_gene.tbl.txt",
  delim ="\t", skip = 29)

# translate feature1 and feature2 ids to TF/target symbols
fantom$tf = fantom_anno$primary_name[match(fantom$feature1_id, fantom_anno$feature_id)]
fantom$target = fantom_anno$primary_name[match(fantom$feature2_id, fantom_anno$feature_id)]


# subset fantom dataset to tfs available in tf census
fantom_clean = fantom %>%
  filter(tf %in% tfs) %>%
  arrange(tf, target) %>%
  distinct(tf, target, pubmed_id) %>%
  mutate(pubmed_id = as.character(pubmed_id))

saveRDS(fantom_clean,
        "inst/extdata/networks/curated/fantom_4/network_with_pubmed.rds")

#### TRRUST ####
trrust = read_delim(
  "inst/extdata/tf_target_sources/curated/trrust/trrust_rawdata.human.v2.tsv",
  delim = "\t", col_names = c("tf", "target", "mor", "pubmed_id")
  )

trrust_clean = trrust %>%
  separate_rows(pubmed_id, sep = ";", ) %>%
  add_count(tf, target, mor, name = "evidence_count") %>%
  group_by(tf, target, mor, evidence_count) %>%
  summarise(pubmed_id = str_c(pubmed_id, collapse = ",")) %>%
  ungroup() %>%
  mutate(priority = case_when(mor == "Unknown" ~ FALSE,
                              TRUE ~ TRUE)) %>%
  arrange(-priority, -evidence_count, mor) %>%
  group_by(tf, target) %>%
  slice(1) %>%
  ungroup() %>%
  arrange(tf, target) %>%
  select(-priority) %>%
  distinct(tf, target, mor, pubmed_id) %>%
  mutate(mor = case_when(mor == "Repression" ~ -1,
                         mor == "Unknown" ~ 0,
                         mor == "Activation" ~ 1)) %>%
  filter(tf %in% tfs)

saveRDS(trrust_clean,
        "inst/extdata/networks/curated/trrust/network_with_pubmed.rds")

#### IntAct ####
int_act = read_delim(
  "inst/extdata/tf_target_sources/curated/int_act/homo_sapiens_protein2gene_20171031.txt",
  delim = "\t"
  ) %>%
  janitor::clean_names()

network = int_act %>%
  select(number_id_s_interactor_a, id_s_interactor_b, type_s_interactor_a,
         publication_identifier_s) %>%
  # assign tf/target class based on the type of interactor a
  mutate(tf = case_when(type_s_interactor_a == "psi-mi:\"MI:0250\"(gene)" ~ id_s_interactor_b,
                        TRUE ~ number_id_s_interactor_a),
         target = case_when(type_s_interactor_a == "psi-mi:\"MI:0250\"(gene)" ~ number_id_s_interactor_a,
                            TRUE ~ id_s_interactor_b)) %>%
  select(tf, target, pubmed_id = publication_identifier_s) %>%
  # separate_rows(pubmed_id, sep = "\\|") %>%
  # filter(str_detect(pubmed_id, "pubmed:")) %>%
  # mutate(pubmed_id = str_remove(pubmed_id, "pubmed:")) %>%
  # ungroup() %>%
  mutate(tf = str_remove(str_remove(tf, "uniprotkb:"), "-[:graph:]*"),
         target = str_remove(target, "ensembl:"))

# map ensembl gene and swissprot ids to HGNC symbols
anno = readRDS("inst/extdata/annotations/gene_annotation.rds")
network$tf = anno$hgnc_symbol[match(network$tf, anno$uniprotswissprot)]
network$target = anno$hgnc_symbol[ match(network$target, anno$ensembl_gene_id)]

# subset tfs to census tfs
int_act_clean = network %>%
  drop_na() %>%
  filter(tf %in% tfs) %>%
  group_by(tf, target) %>%
  summarise(pubmed_id = str_c(unique(pubmed_id), collapse = ",")) %>%
  ungroup() %>%
  arrange(tf, target) %>%
  distinct()

saveRDS(int_act_clean,
        "inst/extdata/networks/curated/int_act/network_with_pubmed.rds")

#### Oreganno ####
df = read_delim(pipe("grep sapiens inst/extdata/tf_target_sources/curated/oreganno/ORegAnno_Combined_2016.01.19.tsv | grep TRANSCRIPTION | grep hg38"),
                delim = "\t", col_names = FALSE)

df_clean = df %>%
  select(source = X13,
         effect = X3,
         Regulatory_Element_Source = X10,
         tf = X8,
         target = X5,
         pubmed_id = X12) %>%
  filter(Regulatory_Element_Source != "miRBase" &
           tf %in% tfs &
           target != "N/A") %>%
  drop_na(target) %>%
  mutate(mor = case_when(effect == "POSITIVE OUTCOME" ~ 1,
                         effect == "NEGATIVE OUTCOME" ~ -1,
                         TRUE ~ 0)) %>%
  select(-effect, -Regulatory_Element_Source) %>%
  arrange(tf)

oreganno = df_clean %>%
  filter(source %in% c("N/A", "STAT1 literature-derived sites")) %>%
  distinct(tf, target, mor, pubmed_id)


regulome = df_clean %>%
  filter(source == "NFIRegulomeDB") %>%
  distinct(tf, target, mor)

pazar = df_clean %>%
  filter(source == "PAZAR") %>%
  distinct(tf, target, mor)

saveRDS(oreganno,
        "inst/extdata/networks/curated/oreganno/network_with_pubmed.rds")
saveRDS(regulome, "inst/extdata/networks/curated/nfi_regulome_db/network.rds")
saveRDS(pazar, "inst/extdata/networks/curated/pazar/network.rds")


#### TFact ####
df = read_delim("inst/extdata/tf_target_sources/curated/tf_act/Catalogues.txt",
                    delim = "\t") %>%
  janitor::clean_names()


df_clean = df %>%
  # filter for human data
  filter(str_detect(species, fixed("homo", ignore_case = TRUE)) |
           str_detect(species, fixed("human", ignore_case = TRUE))) %>%
  select(source = ref, tf = official_tf_coding_gene_name,
         target = official_gene_name, mor = regulation,
         pubmed_id = ref_accession_number) %>%
  mutate(mor = case_when(mor == "UP" ~ 1,
                         mor == "DOWN" ~ -1,
                         TRUE ~ 0)) %>%
  mutate(pubmed_id = str_remove(pubmed_id, pattern = "^;"),
         pubmed_id = str_remove(pubmed_id, pattern = ";$")) %>%
  separate_rows(pubmed_id, sep = ";") %>%
  filter(!str_detect(pubmed_id, "pazar")) %>%
  group_by(source, tf, target, mor) %>%
  summarise(pubmed_id = str_c(pubmed_id, collapse = ",")) %>%
  ungroup()


tf_act = df_clean %>%
  filter(str_detect(source, fixed("pubmed", ignore_case = TRUE))) %>%
  distinct(tf, target, mor, pubmed_id) %>%
  filter(tf %in% tfs) %>%
  arrange(tf, target)

trrd = df_clean %>%
  filter(str_detect(source, fixed("trrd", ignore_case = TRUE))) %>%
  distinct(tf, target, mor) %>%
  filter(tf %in% tfs) %>%
  arrange(tf, target)


saveRDS(tf_act,
        "inst/extdata/networks/curated/tf_act/network_with_pubmed.rds")
saveRDS(trrd, "inst/extdata/networks/curated/trrd_via_tf_act/network.rds")

#### Reviews ####
review_paths = list.files("inst/extdata/tf_target_sources/curated/reviews",
                          pattern = "sif", full.names = TRUE,
                          recursive = TRUE) %>%
  discard(.p = ~str_detect(.x, "network"))

reviews = review_paths %>%
  map_dfr(function(path) {
    pubmed_id = path %>%
      str_split("/") %>%
      pluck(1,6) %>%
      str_split("_") %>%
      pluck(1) %>%
      tail(1)

    read_delim(path, delim = "\t", col_names = FALSE) %>%
      rename(tf = X1, target = X2) %>%
      filter(tf != "TF") %>%
      select(tf, target) %>%
      mutate(pubmed_id = pubmed_id)
  })

saveRDS(reviews,
        "inst/extdata/networks/curated/reviews/network_with_pubmed.rds")

#### TRED ####
tred = read_csv(
  "inst/extdata/tf_target_sources/curated/tred_via_reg_network/export_Fri_Sep_22_11_00_16_UTC_2017.csv")

tred_clean = tred %>%
  distinct(tf = regulator_symbol, target = target_symbol) %>%
  arrange(tf, target) %>%
  filter(tf %in% tfs)

saveRDS(tred_clean,
        "inst/extdata/networks/curated/tred_via_reg_network/network.rds")

#### HTRIdb ####
htri = read_delim("inst/extdata/tf_target_sources/curated/htri_db/tf-target_network_052016_literaturecurated.txt",
                  delim = "\t")

htri_clean = htri %>%
  distinct(tf = SYMBOL_TF, target = SYMBOL_TG, pubmed_id = PUBMED_ID) %>%
  filter(tf %in% tfs) %>%
  arrange(tf, target) %>%
  distinct() %>%
  mutate(pubmed_id = as.character(pubmed_id))

saveRDS(htri_clean,
        "inst/extdata/networks/curated/htri_db/network_with_pubmed.rds")

#### KEGG ####

get_GErel_edges = function(path_name, path_id){

  message(path_name)
  tmp = paste(tempfile(tmpdir = "~/tmp"), ".xml", sep = "")
  my_kgml = retrieveKGML(pathwayid = path_id, organism="hsa", destfile = tmp)
  mapkG = try(parseKGML2Graph(tmp, expandGenes=TRUE), silent = TRUE)

  if (is(mapkG,"try-error"))
    return(NULL)

  pathedges =  getKEGGedgeData(mapkG)
  pathedges_type = sapply(pathedges, getType)
  pathedges_GErel = pathedges[ pathedges_type == "GErel" ]
  pathedges_entryID = sapply(pathedges_GErel, getEntryID)
  pathedges_subtype_name = sapply(pathedges_GErel,
                                  function(x) getName(getSubtype(x)[[1]]))
  pathedges_subtype_value = sapply(pathedges_GErel,
                                   function(x) getValue(getSubtype(x)[[1]]))
  if (is.list(pathedges_entryID))
    return(NULL)
  df = as.data.frame(cbind(t(pathedges_entryID),
                           pathedges_subtype_name,
                           pathedges_subtype_value,
                           path_name = path_name,
                           path_id = path_id),
                     stringsAsFactors = FALSE)

  df$Entry1Symbol = sapply(mget(translateKEGGID2GeneID(df$Entry1ID),
                                org.Hs.egSYMBOL, ifnotfound=NA), "[[",1)
  df$Entry2Symbol = sapply(mget(translateKEGGID2GeneID(df$Entry2ID),
                                org.Hs.egSYMBOL, ifnotfound=NA), "[[",1)
  return(df)
}



mykegglist = keggList(organism = "hsa", database = "pathway")
pId = gsub("path:hsa", "", names(mykegglist))
pName = mykegglist
GErel_edges = mapply(get_GErel_edges, pName, pId)


x = GErel_edges %>%
  compact() %>%
  enframe(name = "key") %>%
  unnest(value) %>%
  arrange(Entry1Symbol, Entry2Symbol)

network = x %>%
  distinct(tf = Entry1Symbol, target = Entry2Symbol) %>%
  filter(tf %in% tfs) %>%
  arrange(tf, target)

saveRDS(network, "inst/extdata/networks/curated/kegg/network.rds")
