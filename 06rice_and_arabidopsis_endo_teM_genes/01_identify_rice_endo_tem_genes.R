######rice expression data in multiple tissues

metadata <- read.csv("~/Documents/lab/endo_tem_genes_paper/other plants/RED_2_0_bulk_metadata.csv",header=T)

rice_expression_raw <- read.delim("~/Documents/lab/endo_tem_genes_paper/other plants/rice_expression_raw.txt", comment.char="#")
tpm <- function(x) (x/rice_expression_raw$Length)/(sum(x/rice_expression_raw$Length))*1e6
rice_tpm <- lapply(rice_expression_raw[7:21], tpm)
rice_tpm<-as.data.frame(rice_tpm)
rice_tpm$gene <- rice_expression_raw$Geneid

colnames(rice_tpm)[1:15]<-sub(
  "^\\.\\.bam\\.(SRR[0-9]+)Aligned\\.sortedByCoord\\.out\\.bam$",
  "\\1",
  colnames(rice_tpm)[1:15]
)

library(dplyr)
library(tidyr)

tpm_matrix<-rice_tpm

subset_metadata<-metadata%>%
  filter(Run%in%colnames(tpm_matrix))%>%
  select(c(Run,tissue_normalized,development_stage))%>%
  mutate(endo=ifelse(tissue_normalized=="endosperm","endo","others"))

tpm_matrix<-tpm_matrix%>%
  rowwise() %>%
  mutate(
    endo_max  = max(c_across(all_of(subset_metadata[which(subset_metadata$endo=="endo"),"Run"])), na.rm = TRUE),
    other_max = max(c_across(all_of(subset_metadata[which(subset_metadata$endo=="others"),"Run"])), na.rm = TRUE),
    fc=(endo_max+1)/(other_max+1)
  ) %>%
  ungroup()


###############################rice leaf methylation data 
rice_methy_nip <- read.csv("~/Documents/lab/endo_tem_genes_paper/other plants/rice_methy_nip.cvs", sep="")

wide_df <- rice_methy_nip %>%
  select(gene, context, mC, umC,totalC, methy) %>%
  pivot_wider(
    names_from = context,
    values_from = c(mC, umC,totalC, methy),
    names_glue = "{context}_{.value}"
  )


rice.methylation_leaf <- wide_df %>%
  mutate(epiallele = case_when(
    CG_totalC==0 & CHG_totalC ==0 & CHH_totalC==0 ~ "no_coveraga",
    CG_methy <= 0.05 & CHG_methy <= 0.05  ~  "UM",
    CG_methy >= 0.2 & CHG_methy <= 0.05  ~  "gbM",
    CG_methy >= 0.4 & CHG_methy >= 0.4 & CHG_totalC >= 10 ~  "teM",
    TRUE ~ "ambiguous"         
  ))


rice_methy_tpm<-merge(rice.methylation_leaf,tpm_matrix,by="gene")%>%
  select(c(1:14,30:32))
rice_alltem_genes<-rice_methy_tpm%>%
  filter(epiallele=="teM"&endo_max>=20)

rice_expressed_genes<-rice_methy_tpm %>%
  filter(endo_max>=20)

rice_specific_genes<-rice_methy_tpm %>%
  filter(endo_max >=20 & fc >=5)

rice_tem_specific_genes<-rice_specific_genes%>%
  filter(epiallele=="teM")

