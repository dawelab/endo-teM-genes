#identify imprinted genes using the muli-line data
mat_pref_NAM <- read.delim("~/Documents/lab/megpeg/mat_pref_NAM.txt")
MaizeGDB_maize_pangene_2020_08 <- read.delim("~/Documents/lab/MaizeGDB_maize_pangene_2020_08.tsv", header=FALSE)
library(tidyr)
library(dplyr)

pangene_long <- MaizeGDB_maize_pangene_2020_08 %>%
  mutate(pangene_id = row_number()) %>%
  pivot_longer(-pangene_id, names_to = "col", values_to = "feature") %>%
  filter(!is.na(feature))

pangene_b73_map <- pangene_long %>%
  group_by(pangene_id) %>%
  mutate(b73_gene = first(feature[grepl("^Zm00001eb", feature)])) %>%
  ungroup() %>%
  select(feature, b73_gene)%>%
  filter(!is.na(feature))

pangene_b73_map <- pangene_long %>%
  group_by(pangene_id) %>%
  filter(any(grepl("^Zm00001eb", feature))) %>%
  summarise(
    b73_gene = feature[grepl("^Zm00001eb", feature)],
    gene = list(feature),
    .groups = "drop"
  ) %>%
  unnest(gene) %>%
  select(feature = gene, b73_gene)

mat_pref <- mat_pref_NAM %>%
  left_join(pangene_b73_map, by = "feature") %>%
  filter(contrast != "OB_MN") %>%
  mutate(
    maternal_preference=case_when(A_mean<=5 & B_mean <=5~ NA_real_, TRUE ~ maternal_preference ))
  
mat_pref_var<-mat_pref%>%
  group_by(b73_gene)%>%
  mutate(variance=var(maternal_preference,na.rm = TRUE))%>%
  ungroup()%>%
  filter(variance<=0.02)
