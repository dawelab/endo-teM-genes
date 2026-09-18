ara.expression_endo_tem<-endosperm_specific_expressed%>%
  filter(gene%in%ara_endo_tem_genes)%>%
  mutate(class="endo teM genes")%>%
  select(c("gene","endo_max","fc","class"))

ara.expression_endo_specific <- endosperm_specific_expressed%>%
  mutate(class="endosperm genes")%>%
  select(c("gene","endo_max","fc","class"))
ara.expression_endo_expressed<-endosperm_highly_expressed%>%
  mutate(class="endo expressed genes")%>%
  select(c("gene","endo_max","fc","class"))
ara.expression.groups<-rbind( ara.expression_endo_specific,ara.expression_endo_tem,ara.expression_endo_expressed)

ara.expression.groups$class <- factor(
  ara.expression.groups$class,
  levels = c(
    "endo teM genes",
    "endosperm genes",
    "endo expressed genes"
  )
)
representative_transcript <- read.delim("~/Documents/lab/endo_tem_genes_paper/other plants/arabidopsis/representative_transcript.txt", header=FALSE)
representative_transcript$id<-substring(representative_transcript$V1,2,10)
ara.expression.groups<-ara.expression.groups%>%
  filter(gene%in%representative_transcript$id)
ggplot(
  ara.expression.groups,
  aes(x = log2(fc), y = log2(endo_max))
) +
  geom_point(
    data = subset(ara.expression.groups, class == "endo expressed genes"),
    aes(color = class),
    alpha = 0.9,
    size = 1.2
  ) +
  geom_point(
    data = subset(ara.expression.groups, class == "endosperm genes"),
    aes(color = class),
    alpha = 0.9,
    size = 1.2
  ) +
  geom_point(
    data = subset(ara.expression.groups, class == "endo teM genes"),
    aes(color = class),
    alpha = 0.9,
    size = 1.2
  ) +
  scale_color_manual(
    values = c(
      "endo expressed genes" = "gray80",
      "endosperm genes" = "salmon",
      "endo teM genes" = "deepskyblue"
    ),
    breaks = c(
      "endo teM genes",
      "endosperm genes",
      "endo expressed genes"
    )
  ) +
  theme_classic(base_size = 16) +
  labs(
    x = "log2(FC)",
    y = "log2(CPM of endosperm)",
    color = NULL
  ) 


######### bar plots of methylation in leaf, wt endo, dme mutant endo
merged_all_methy <- read.csv("~/Documents/lab/endo_tem_genes_paper/other plants/arabidopsis/dme methylation/merged_all_methy.cvs", sep="")

library(dplyr)
library(tidyr)
ara.methylation_endo <- merged_all_methy %>%
  select(gene, sample, context, mC, umC, totalC, methy) %>%
  pivot_wider(
    #id_cols = c(gene, sample),
    names_from = context,
    values_from = c(mC, umC, totalC, methy),
    names_glue = "{context}_{.value}"
  )

ara.methylation_leaf$sample<-"leaf"
ara.methylation_leaf2<-ara.methylation_leaf%>%
  select(1,15,2:13)

ara.methylation_leaf_endo<-rbind(ara.methylation_leaf2,ara.methylation_endo)

ara_endo_tem_genes<-endosperm_specific_expressed$gene[which(endosperm_specific_expressed$epiallele=="teM")]

ara.methylation_endo_tem<-ara.methylation_leaf_endo%>%
  filter(gene%in%ara_endo_tem_genes)%>%
  mutate(class="endo teM genes")

ara.methylation_endo_specific<-ara.methylation_leaf_endo%>%
  filter(gene%in%endosperm_specific_expressed$gene)%>%
  mutate(class="endosperm specific genes")

ara.methylation_endo_expressed<-ara.methylation_leaf_endo%>%
  filter(gene%in%endosperm_highly_expressed$gene)%>%
  mutate(class="endosperm expressed genes")

ara.methylation.groups<-rbind(ara.methylation_endo_tem,ara.methylation_endo_specific,ara.methylation_endo_expressed)


library(ggplot2)

ara.methylation.groups$class <- factor(
  ara.methylation.groups$class,
  levels = c("endo teM genes", "endosperm specific genes", "endosperm expressed genes")
)

ara.methylation.groups$sample <- factor(
  ara.methylation.groups$sample,
  levels = c("leaf", "wt", "dme")
)

ggplot(ara.methylation.groups, aes(
  x = class,
  y = CG_methy,
  fill = sample
)) +
  geom_boxplot(
    position = position_dodge(width = 0.7),
    width = 0.5,
    alpha = 0.7,
    outliers = T
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
      size = 15
    ),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.title = element_blank(),
    legend.text = element_text(size = 15)
  )


mCG_wide <- ara.methylation.groups %>%
  filter(
    sample %in% c("wt", "leaf"),
    !is.na(CG_methy)
  ) %>%
  select(gene, class, sample, CG_methy) %>%
  distinct() %>%
  pivot_wider(
    names_from = sample,
    values_from = CG_methy
  ) %>%
  filter(
    !is.na(wt),
    !is.na(leaf)
  )

mCG_test <- mCG_wide %>%
  group_by(class) %>%
  summarise(
    n_pairs = n(),
    endosperm_median = median(wt),
    leaf_median = median(leaf),
    median_difference = median(wt - leaf),
    p_value = wilcox.test(
      wt,
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

mCG_wide <- ara.methylation.groups %>%
  filter(
    sample %in% c("wt", "dme"),
    !is.na(CG_methy)
  ) %>%
  select(gene, class, sample, CG_methy) %>%
  distinct() %>%
  pivot_wider(
    names_from = sample,
    values_from = CG_methy
  ) %>%
  filter(
    !is.na(wt),
    !is.na(dme)
  )

mCG_test <- mCG_wide %>%
  group_by(class) %>%
  summarise(
    n_pairs = n(),
    wt_median = median(wt),
    dme_median = median(dme),
    median_difference = median(dme - wt),
    p_value = wilcox.test(
      dme,
      wt,
      paired = TRUE,
      alternative = "greater",
      exact = FALSE
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_adjusted = p.adjust(p_value, method = "BH")
  )

mCG_test



######ara imprinting
ara_imprinting <- read_csv("Documents/lab/endo_tem_genes_paper/other plants/arabidopsis/ara_imprinting.csv")


ara_imprinting <-ara_imprinting %>%
  count(`Gene ID`, `Overall Imprinting Status`) %>%
  pivot_wider(
    names_from = `Overall Imprinting Status`,
    values_from = n,
    values_fill = 0
  )%>%
  mutate(
    imprinting_overall = case_when(
      MEG >= 3 & PEG == 0 ~ "MEG",
      PEG >= 3 & MEG == 0 ~ "PEG",
      `not imprinted` >= 3 ~ "not imprinted",
      TRUE ~ "no data"
    ))


ara_expresed_imprinting<-merge(ara_imprinting,endosperm_highly_expressed,by.x="Gene ID",by.y="gene",all.y=T)
ara_specific_imprinting<-merge(ara_imprinting,endosperm_specific_expressed,by.x="Gene ID",by.y="gene",all.y=T)

endo_tem<-endosperm_specific_expressed[which(endosperm_specific_expressed$epiallele=="teM"),]
ara_endo_tem_imprinting<-merge(ara_imprinting,endo_tem,by.x="Gene ID",by.y="gene",all.y=T)






