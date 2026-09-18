#arabidopsis expression data from single cell in seeds and bulk-rna seq in other tissues


library(Seurat)

obj <- readRDS("~/Documents/lab/endo_tem_genes_paper/GSE295007_ATLAS_merged_annotated_sigmods.rds")


raw_counts <- GetAssayData(obj, assay = "RNA", layer = "counts")


groups <- unique(obj$level_1_annotation_timed)

pseudo_bulk_counts <- sapply(groups, function(g) {
  cells <- colnames(obj)[obj$level_1_annotation_timed == g]
  Matrix::rowSums(raw_counts[, cells, drop = FALSE])
})

pseudo_bulk_counts <- as.matrix(pseudo_bulk_counts)

pseudobulk_CPM_matrix<- sweep(
  pseudo_bulk_counts,
  2,
  colSums(pseudo_bulk_counts),
  FUN = "/"
) * 1e6





library(readr)
E_MTAB_7978_query_results_tpmss <- read_delim("Documents/lab/endo_tem_genes_paper/other plants/arabidopsis/E-MTAB-7978-query-results.tpmss.tsv", 
                                              delim = "\t", escape_double = FALSE, 
                                              comment = "#", trim_ws = TRUE)


all_tissue_tpm <- E_MTAB_7978_query_results_tpmss %>%
  select(-c("dry seed stage, seed, not applicable",
            "seed imbibition stage, seed, not applicable",
            "silique stage 1, silique, not applicable",
            "silique stage 2, silique, not applicable",
            "silique stage 3, silique, not applicable",
            "silique stage 4, silique, not applicable",
            "silique stage 5, silique, not applicable",
            "callus, plant callus, not applicable" )) %>%
  rowwise() %>%
  mutate(
    atlas_other_max = max(c_across(3:48), na.rm = TRUE),
  ) %>%
  ungroup()

all_tissue_tpm  <- all_tissue_tpm  %>%
  mutate(
    atlas_other_max  = replace(atlas_other_max , is.infinite(atlas_other_max), NA)
  )
library(tibble)

pseudobulk_CPM_matrix <- pseudobulk_CPM_matrix %>%
  as.data.frame() %>%
  rownames_to_column("gene")

endosperm_highly_expressed<-pseudobulk_CPM_matrix %>%
  rowwise() %>%
  mutate(endo_max=max(c_across(c("3DAP Endosperm","5DAP Endosperm","7DAP Endosperm"))),
         other_tissue_max=max(c_across(-c("gene","3DAP Endosperm","5DAP Endosperm","7DAP Endosperm","endo_max")))) %>%
  ungroup()%>%
  mutate(fc=(endo_max+1)/(other_tissue_max+1)) %>%
  filter(endo_max >= 20)

####methylation data
ara_methy <- read.csv("~/Documents/lab/endo_tem_genes_paper/other plants/arabidopsis/ara_methy.cvs", sep="")

library(dplyr)
library(tidyr)
wide_df <- ara_methy %>%
  select(gene, context, mC, umC,totalC, methy) %>%
  pivot_wider(
    names_from = context,
    values_from = c(mC, umC,totalC, methy),
    names_glue = "{context}_{.value}"
  )



ara.methylation_leaf <- wide_df %>%
  mutate(epiallele = case_when(
    CG_totalC==0 & CHG_totalC ==0 & CHH_totalC==0 ~ "no_coveraga",
    CG_methy <= 0.05 & CHG_methy <= 0.05  ~  "UM",
    CG_methy >= 0.2 & CHG_methy <= 0.05  ~  "gbM",
    CG_methy >= 0.4 & CHG_methy >= 0.2 & CHG_totalC >= 10 ~  "teM",
    TRUE ~ "ambiguous"         
  ))



endosperm_highly_expressed<-merge(endosperm_highly_expressed,ara.methylation_leaf,by="gene",all.x=T,al.y=F)
ara_alltem<-endosperm_highly_expressed%>%
  filter(epiallele=="teM")

endosperm_specific_expressed0<-merge(endosperm_highly_expressed,all_tissue_tpm,by.x="gene",by.y="Gene ID")%>%
  filter(fc>=5 )
endosperm_specific_expressed<-endosperm_specific_expressed0%>%
  filter(atlas_other_max<=4)



ara_endo_tem_genes<-endosperm_specific_expressed$gene[which(endosperm_specific_expressed$epiallele=="teM")]
