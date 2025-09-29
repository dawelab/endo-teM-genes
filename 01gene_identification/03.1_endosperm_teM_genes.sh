#identify teM genes, mCG >= 0.4 & mCHG >= 0.4 & cCHG >= 10 
B73.all <- read.csv("~/Documents/lab/megpeg/B73.all.csv")
all_tpm <- read.csv("~/Documents/lab/all_tpm.txt", sep="")

b73.all.6_38<-merge(B73.all,all_tpm,by="gene",all=T)
b73.all.6_38 <- b73.all.6_38 %>%
  mutate(epiallele2 = case_when(
    mCG <= 0.05 & mCHG <= 0.05  ~  "UM",
    mCG >= 0.2 & mCHG <= 0.05  ~  "gbM",
    mCG >= 0.4 & mCHG >= 0.4 & cCHG >= 10 ~  "teM",
    TRUE ~ "ambiguous"         
  ))
  
#identify endoseprm teM genes,TPM >= 5fold and tpm >= 20 in endosperm#############

all_body_genes<-data.frame(matrix(nrow = 0, ncol = 1))

for (i in seq(6, 38, by = 2)) {
  sample<-paste("D",i,sep = "")
  selected <- b73.all.6_38[apply(b73.all.6_38, 1, function(selected_row) {
    as.numeric(selected_row[[sample]]) >= 5 *max(as.numeric(selected_row[c(14,16:23)]))
  }), ]
  #identify the teM genes for each development stage
  b73.one_stage <- selected %>%
    mutate(zeinlike = case_when(
      epiallele2 == "teM"  & !!sym(sample) >= 20 ~  "TE-like-genes",
      TRUE ~ "none"         
    ))
  
  body_genes<-b73.one_stage[c(which(b73.one_stage$zeinlike=="TE-like-genes")),1]
  
  # name<-paste("~/Documents/lab/megpeg/BM_5FC",i, ".bed", sep = "")
  # write.table(body_genes,name,quote = F,row.names = F,
  #             col.names = F,sep = "\t")
  
  #combind all the stages
  all_body_genes<-rbind(all_body_genes,as.data.frame(body_genes))
  
}

all_body_genes<-as.data.frame(unique(all_body_genes$body_genes))
colnames(all_body_genes)<-"teM_genes"
# name<-paste("~/Documents/lab/megpeg/all_BM_5FC.bed", sep = "")
# write.table(all_body_genes,name,quote = F,row.names = F,
#            col.names = F,sep = "\t")

colnames(all_body_genes)<-'b73'

#core genes
b73_core_genes <- read.csv("~/Documents/lab/b73_core_genes.txt", sep="")
core_gene_te_like_gene<-merge(b73_core_genes,all_body_genes,by.y = "b73",by.x = "x", all=F)
core_gene_te_like_gene<-core_gene_te_like_gene$x
#write.table(core_gene_te_like_gene,"~/Documents/lab/megpeg/endosperm_teM_genes.txt",quote = F,row.names = F,
#            col.names = F,sep = "\t")
