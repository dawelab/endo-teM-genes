#TPM >= 5fold and tpm >= 20 in endosperm#############

all_body_genes<-data.frame(matrix(nrow = 0, ncol = 1))

for (i in seq(6, 38, by = 2)) {
  sample<-paste("D",i,sep = "")
  selected <- b73.all.6_38[apply(b73.all.6_38, 1, function(selected_row) {
    as.numeric(selected_row[[sample]]) >= 5 *max(as.numeric(selected_row[c(14,16:23)]))
  }), ]
  
  b73.one_stage <- selected %>%
    mutate(zeinlike = case_when(
      epiallele2 == "teM"  & !!sym(sample) >= 20 ~  "TE-like-genes",
      TRUE ~ "none"         
    ))
  
  body_genes<-b73.one_stage[c(which(b73.one_stage$zeinlike=="TE-like-genes")),1]
  
  # name<-paste("~/Documents/lab/megpeg/BM_5FC",i, ".bed", sep = "")
  # write.table(body_genes,name,quote = F,row.names = F,
  #             col.names = F,sep = "\t")
  
  all_body_genes<-rbind(all_body_genes,as.data.frame(body_genes))
  
}

all_body_genes<-as.data.frame(unique(all_body_genes$body_genes))
colnames(all_body_genes)<-"teM_genes"
# name<-paste("~/Documents/lab/megpeg/all_BM_5FC.bed", sep = "")
# write.table(all_body_genes,name,quote = F,row.names = F,
#            col.names = F,sep = "\t")

colnames(all_body_genes)<-'b73'
core_gene_te_like_gene<-merge(core_gene,all_body_genes,by.y = "b73",by.x = "gene", all=F)
core_gene_te_like_gene<-core_gene_te_like_gene$gene
#write.table(core_gene_te_like_gene,"~/Documents/lab/megpeg/all_BM_5FC_core.bed",quote = F,row.names = F,
#            col.names = F,sep = "\t")

