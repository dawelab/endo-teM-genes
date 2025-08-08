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


bmegs_var<-mat_pref_var%>%
  filter(b73_gene%in%all_BM_5FC_core$V1)%>%
  group_by(b73_gene)%>%
  summarise(mean=mean(maternal_preference,na.rm=T))%>%
  mutate(group="BMEG")


fmegs_var<-mat_pref_var%>%
  filter(b73_gene%in%MEG_FMEGs$V1)%>%
  group_by(b73_gene)%>%
  summarise(mean=mean(maternal_preference,na.rm=T))%>%
  mutate(group="FMEG")

coregene_var<-mat_pref_var%>%
  filter(b73_gene%in%core_gene$gene)%>%
  group_by(b73_gene)%>%
  summarise(mean=mean(maternal_preference,na.rm=T))%>%
  mutate(group="Core Gene")

combined_var<-rbind(bmegs_var,fmegs_var,coregene_var)%>%
  mutate(group=factor(group,levels=c("BMEG","FMEG","Core Gene")))

ggplot(combined_var, aes(x = group, y = mean,color=group)) +
  geom_point() +
  #geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  theme_minimal() +
  theme(text = element_text(size = 18))+
  labs(title = "Maternal expression  preference ",
       x = "",
       y = "Maternal Expression Preference") +
  geom_hline(yintercept = 0.86, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0.4, linetype = "dashed", color = "black")
