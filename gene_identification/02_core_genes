pan_gene_matrix_v3_cyverse <- read.csv("~/Desktop/pan_gene_matrix_v3_cyverse.csv")

library(dplyr)
library(tidyr)
library(stringr)

#separate the B73 genes corrosponding to the same pangene id
df_long <- pan_gene_matrix_v3_cyverse %>%
  separate_rows(B73, sep = ";\\s*") %>%
  mutate(b73_gene = substring(B73,1,15))

core_genes<-df_long%>%
  filter(class=="Core Gene")%>%
  select(b73_gene,class)

#request the genes are identified in B73 v5 reference genome
Zm.B73.REFERENCE<- read.delim("~/Documents/lab/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3", comment.char="#")
b73_genes<-Zm.B73.REFERENCE%>%
  filter(chromosome=='gene')%>%
  mutate(gene=substring(ID.1.Name.chromosome.Zm.B73.REFERENCE.NAM.5.0.1.1.308452471.1,4,18))

b73_core_genes<-merge(core_genes,b73_genes,by.x="b73_gene",by.y="gene")

write.csv(b73_core_genes$b73_gene,"~/Documents/lab/b73_core_genes.txt",quote=F,row.names = F,col.names = F)  

