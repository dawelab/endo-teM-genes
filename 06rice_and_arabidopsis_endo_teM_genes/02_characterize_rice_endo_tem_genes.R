#########boxplots of methylation in leaf and endosperm

rice_methy_leaf_endo <- read.csv("~/Documents/lab/endo_tem_genes_paper/other plants/rice_methy_leaf_endo.cvs", sep="") # nolint # nolint: line_length_linter.

wide_df <-rice_methy_leaf_endo %>%
  select(gene, context, mC, umC,totalC, methy,sample) %>%
  pivot_wider(
    names_from = context,
    values_from = c(mC, umC,totalC, methy),
    names_glue = "{context}_{.value}"
  )

rice.methylation_endo <- wide_df  %>% # nolint
  filter(sample=="SRR9637047")%>%
  mutate(epiallele = case_when(
    CG_totalC==0 & CHG_totalC ==0 & CHH_totalC==0 ~ "no_coveraga",
    CG_methy <= 0.05 & CHG_methy <= 0.05  ~  "UM",
    CG_methy >= 0.2 & CHG_methy <= 0.05  ~  "gbM",
    CG_methy >= 0.4 & CHG_methy >= 0.4 & CHG_totalC >= 10 ~  "teM",
    TRUE ~ "ambiguous"         
  ))%>%
  mutate(tissue="endosperm")

rice.methylation_endo_tem<-rice.methylation_endo%>%
  filter(gene%in%rice_tem_specific_genes$gene)%>%
  mutate(class="endo teM genes")

rice.methylation_endo_specific<-rice.methylation_endo%>%
  filter(gene%in%rice_specific_genes$gene)%>%
  mutate(class="endosperm specific genes")

rice.methylation_endo_expressed<-rice.methylation_endo%>%
  filter(gene%in%rice_expressed_genes$gene)%>%
  mutate(class="endosperm expressed genes")

rice.methylation.groups<-rbind(rice.methylation_endo_tem,rice.methylation_endo_specific,rice.methylation_endo_expressed)



rice.methylation_leaf <- wide_df %>%
  filter(sample=="SRR9637033")%>%
  mutate(epiallele = case_when(
    CG_totalC==0 & CHG_totalC ==0 & CHH_totalC==0 ~ "no_coveraga",
    CG_methy <= 0.05 & CHG_methy <= 0.05  ~  "UM",
    CG_methy >= 0.2 & CHG_methy <= 0.05  ~  "gbM",
    CG_methy >= 0.4 & CHG_methy >= 0.4 & CHG_totalC >= 10 ~  "teM",
    TRUE ~ "ambiguous"         
  ))%>%
  mutate(tissue="leaf")

rice.methylation_leaf_endo_tem<-rice.methylation_leaf%>%
  filter(gene%in%rice_tem_specific_genes$gene)%>%
  mutate(class="endo teM genes")

rice.methylation_leaf_specific<-rice.methylation_leaf%>%
  filter(gene%in%rice_specific_genes$gene)%>%
  mutate(class="endosperm specific genes")

rice.methylation_leaf_expressed<-rice.methylation_leaf%>%
  filter(gene%in%rice_expressed_genes$gene)%>%
  mutate(class="endosperm expressed genes")

rice.methylation.groups<-rbind(rice.methylation.groups,rice.methylation_leaf_endo_tem,rice.methylation_leaf_specific,rice.methylation_leaf_expressed)


library(ggplot2)

rice.methylation.groups$class <- factor(
  rice.methylation.groups$class,
  levels = c("endo teM genes", "endosperm specific genes", "endosperm expressed genes")
)

rice.methylation.groups$tissue<-factor(
  rice.methylation.groups$tissue,
  levels = c("leaf","endosperm")
)
ggplot(rice.methylation.groups, aes(
  x = class,
  y = CG_methy,
  fill = tissue
)) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.65,
    alpha = 0.7,
    outliers = TRUE
  ) +
  labs(
    x = NULL,
    y = "mCG"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 30,
      hjust = 1,
      size = 12
    ),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  )

library(dplyr)
library(tidyr)
library(rstatix)


mCG_wide <- rice.methylation.groups %>%
  filter(
    tissue %in% c("endosperm", "leaf"),
    !is.na(CG_methy)
  ) %>%
  select(gene, class, tissue, CG_methy) %>%
  distinct() %>%
  pivot_wider(
    names_from = tissue,
    values_from = CG_methy
  ) %>%
  filter(
    !is.na(endosperm),
    !is.na(leaf)
  )

mCG_test <- mCG_wide %>%
  group_by(class) %>%
  summarise(
    n_pairs = n(),
    endosperm_median = median(endosperm),
    leaf_median = median(leaf),
    median_difference = median(endosperm - leaf),
    p_value = wilcox.test(
      endosperm,
      leaf,
      paired = TRUE,
      alternative = "less",
      exact = FALSE
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_adjusted = p.adjust(p_value, method = "BH")
  )

mCG_test


#dot plots of different gene groups

rice.expression_tem<-rice_tem_specific_genes%>%
  mutate(class="endo teM genes")%>%
  select(c("gene","endo_max","fc","class"))


rice.expression_endo_specific <- rice_specific_genes%>%
  mutate(class="endosperm genes")%>%
  select(c("gene","endo_max","fc","class"))
rice.expression_endo_expressed<-rice_expressed_genes%>%
  mutate(class="endo expressed genes")%>%
  select(c("gene","endo_max","fc","class"))

rice.expression_protein_storage<-rice_protein_storage_genes%>%
  select(c(1,4,6))%>%
  mutate(class="protein storage genes")
colnames(rice.expression_protein_storage)[1]<-"gene"
rice.expression.groups<-rbind( rice.expression_protein_storage,rice.expression_endo_specific,rice.expression_tem,rice.expression_endo_expressed)

rice.expression.groups$class <- factor(
  rice.expression.groups$class,
  levels = c(
    "endo teM genes",
    "endosperm genes",
    "endo expressed genes",
    "protein storage genes"
  )
)

ggplot(
  rice.expression.groups,
  aes(x = log2(fc), y = log2(endo_max ))
) +
  geom_point(
    data = subset(rice.expression.groups, class == "endo expressed genes"),
    aes(color = class),
    alpha = 0.9,
    size = 1.2
  ) +
  geom_point(
    data = subset(rice.expression.groups, class == "endosperm genes"),
    aes(color = class),
    alpha = 0.9,
    size = 1.2
  ) +
  geom_point(
    data = subset(rice.expression.groups, class == "endo teM genes"),
    aes(color = class),
    alpha = 0.9,
    size = 1.2
  ) +
  geom_point(
    data = subset(rice.expression.groups, class == "protein storage genes"),
    aes(color = "protein storage genes"),
    shape = 1,
    size = 1.5,
    stroke = 1
  ) +
  scale_color_manual(
    values = c(
      "endo expressed genes" = "gray80",
      "endosperm genes" = "salmon",
      "endo teM genes" = "deepskyblue",
      "protein storage genes" = "black"
    ),
    breaks = c(
      "endo teM genes",
      "endosperm genes",
      "endo expressed genes",
      "protein storage genes"
    )
  ) +
  theme_classic(base_size = 16) +
  labs(
    x = "log2(FC)",
    y = "log2(TPM of endosperm)",
    color = NULL
  )

rice_megs <- read_excel("Documents/lab/endo_tem_genes_paper/other plants/rice imprint classes for MEGs and PEGs..xlsx", 
                        sheet = "MEGs")

rice_pegs <- read_excel("Documents/lab/endo_tem_genes_paper/other plants/rice imprint classes for MEGs and PEGs..xlsx", 
                        sheet = "PEGs")
rice_megs <- rice_megs %>%
  mutate(status = "MEG") %>%
  filter(
    `Imprint class` %in% c(
      "3 and 5 DAP overlapped",
      "3 and 7 DAP overlapped",
      "5 and 7 DAP overlapped",
      "Persistent"
    )
  )
rice_pegs<-rice_pegs%>%
  mutate(status="PEG")%>%
  filter(
    `Imprint class` %in% c(
      "3 and 5 DAP overlapped",
      "3 and 7 DAP overlapped",
      "5 and 7 DAP overlapped",
      "Persistent"
    )
  )

rice_imprinting_info<-rbind(rice_megs,rice_pegs) 

rice_tem_imprinting<-merge(rice_tem_specific_genes,rice_imprinting_info,by.x="gene",by.y="locus name",all.x=T)

rice_expressed_imprinting<-merge(rice_expressed_genes,rice_imprinting_info,by.x="gene",by.y="locus name",all.x=T)

rice_specific_imprinting<-merge(rice_specific_genes,rice_imprinting_info,by.x="gene",by.y="locus name",all.x=T)
