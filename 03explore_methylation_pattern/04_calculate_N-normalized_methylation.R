library(Biostrings)
library(dplyr)

count_methylation_gene <- function(name, sequence) {
  sequence <- toupper(sequence)
  seq_rc <- as.character(reverseComplement(DNAString(sequence)))
  gene_start <- 3
  gene_end   <- nchar(sequence) - 2
  count_CG <- function(seq) {
    sum(sapply(gene_start:gene_end, function(i) {
      if (i <= nchar(seq) - 1) {
        substr(seq, i, i) == "C" && substr(seq, i + 1, i + 1) == "G"
      } else FALSE
    }))
  }
  
  count_CHG <- function(seq) {
    sum(sapply(gene_start:gene_end, function(i) {
      if (i <= nchar(seq) - 2) {
        b1 <- substr(seq, i, i)
        b2 <- substr(seq, i + 1, i + 1)
        b3 <- substr(seq, i + 2, i + 2)
        b1 == "C" && b2 %in% c("A","T","C") && b3 == "G"
      } else FALSE
    }))
  }
  
  count_CHH <- function(seq) {
    sum(sapply(gene_start:gene_end, function(i) {
      if (i <= nchar(seq) - 2) {
        b1 <- substr(seq, i, i)
        b2 <- substr(seq, i + 1, i + 1)
        b3 <- substr(seq, i + 2, i + 2)
        b1 == "C" && b2 %in% c("A","T","C") && b3 %in% c("A","T","C")
      } else FALSE
    }))
  }
  
  
  CG_total  <- count_CG(sequence)  + count_CG(seq_rc)
  CHG_total <- count_CHG(sequence) + count_CHG(seq_rc)
  CHH_total <- count_CHH(sequence) + count_CHH(seq_rc)
  
  tibble(
    gene = name,
    length_bp = nchar(sequence)-4,
    CG_count  = CG_total,
    CHG_count = CHG_total,
    CHH_count = CHH_total
  )
}


#count cytosines in 200bp regions centered TSS
TSS_200bp <- readDNAStringSet("./TSS_200bp.fa", format = "fasta")
seq_names <- names(TSS_200bp)
seqs <- as.character(TSS_200bp)
all_results_TSS_200bp<- purrr::map2_dfr(seq_names, seqs, count_methylation_gene )


write.table(all_results_TSS_200bp"./tss_counts_doublestrand_200bp.csv",quote = F,
            col.names = T,row.names = F,sep = "\t")


           
#count cytosines in gene body regions
gene_body <- readDNAStringSet("./gene_body.fa", format = "fasta")
seq_names <- names(gene_body)
seqs <- as.character(gene_body)
all_results_gene_body<- purrr::map2_dfr(seq_names, seqs, count_methylation_gene )


write.table(all_results_gene_body,"./counts_doublestrand_genebody.csv",quote = F,
            col.names = T,row.names = F,sep = "\t")




#Calculate N-normalized methylation at 200bp centered TSS

#read the C-normalized methylation data
B73_NAM.CGmap.gz.core_gene_B73_10sort_up200.CG.mtr <- read.delim("~/Documents/lab/TSS_methylation/up200_B73/B73_B73v5_merge.CGmap.gz.core_gene_B73_10sort_up200.CG.mtr.gz", header=FALSE)
colnames(B73_NAM.CGmap.gz.core_gene_B73_10sort_up200.CG.mtr)[4]<-"CG"
B73_NAM.CGmap.gz.core_gene_B73_10sort_up200.CHG.mtr <- read.delim("~/Documents/lab/TSS_methylation/up200_B73/B73_B73v5_merge.CGmap.gz.core_gene_B73_10sort_up200.CHG.mtr.gz", header=FALSE)
colnames(B73_NAM.CGmap.gz.core_gene_B73_10sort_up200.CHG.mtr)[4]<-"CHG"

#replace the NA with 0
b73_methylation<-cbind(B73_NAM.CGmap.gz.core_gene_B73_10sort_up200.CG.mtr,B73_NAM.CGmap.gz.core_gene_B73_10sort_up200.CHG.mtr)%>%
  select(c(1,2,3,4,5,11,12))%>%
  mutate(
    CG  = replace_na(as.numeric(CG), 0),
    CHG = replace_na(as.numeric(CHG), 0)
  )

#match the gene ID based on the TSS location
Zm.B73.REFERENCE <- read.delim("~/Documents/lab/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3", header=FALSE, comment.char="#")

TSS_start_site<-Zm.B73.REFERENCE%>%
  filter(V3=="gene")%>%
  mutate(TSS_start_site=ifelse(V7=="+",V4-100-1,V5-100))%>%
  mutate(TSS_start_site=ifelse(TSS_start_site<0,0,TSS_start_site))%>%
  select(1,4,5,7,9,10)%>%
  mutate(V1=substring(V1,4,))

merged <- merge(
  b73_methylation,
  TSS_start_site,
  by.x = c("V1", "V2"),
  by.y = c("V1", "TSS_start_site"),
  all.x = TRUE  )%>%
  mutate(gene=substring(V9,4,18))%>%
  select(gene,CG,CHG)

#merge with the number of cytosine data
tss_counts_doublestrand_200bp <- read.delim("~/Documents/lab/megpeg/tss_counts_doublestrand_200bp.csv")
tss_counts_doublestrand_200bp<-tss_counts_doublestrand_200bp%>%
  mutate(gene=substring(gene, 1, 15))

#caculate the N-normalized methylation
weighted_methylation<-merge(merged, tss_counts_doublestrand_200bp,by="gene")%>%
  mutate(weighted_CG = CG*CG_count/402, # double strand of 200bp centered TSS include the TSS itself.
         weighted_CHG = CHG*CHG_count/402) 

#pollen tem genes
pollen_teM_genes <- read_excel("Documents/lab/TSS_methylation/up200_B73/pollen teM genes.xlsx", 
                               sheet = "candidate DNG target genes")

pollen_teM_genes_nm<-merge(pollen_teM_genes,weighted_methylation,by.x="gene ID",by.y="gene")%>%
  mutate(class=ifelse(is.na(`DEG?`), "non_DEG", "DEG"))%>%
  select("gene ID","weighted_CG","weighted_CHG","class")

colnames(pollen_teM_genes_nm)[1]<-"gene"

#mean(pollen_teM_genes_nm$weighted_CG[which(pollen_teM_genes_nm$class=="DEG")])
#mean(pollen_teM_genes_nm$weighted_CG[which(pollen_teM_genes_nm$class=="non_DEG")])

#mean(pollen_teM_genes_nm$weighted_CHG[which(pollen_teM_genes_nm$class=="DEG")])
#mean(pollen_teM_genes_nm$weighted_CHG[which(pollen_teM_genes_nm$class=="non_DEG")])



#endosperm teM genes 


imprinted_endosperm_tem_genes<-merge(weighted_methylation,bmegs_var,by.x="gene",by.y="b73_gene")

imprinted_endosperm_tem_genes<-imprinted_endosperm_tem_genes%>%
  mutate(class=ifelse(mean>0.85, "maternal", "biparental"))%>%
  select("gene","weighted_CG","weighted_CHG","class")

colnames(imprinted_endosperm_tem_genes)[1]<-"gene"

#mean(imprinted_endosperm_tem_genes$weighted_CG[which(imprinted_endosperm_tem_genes$class=="maternal")])
#mean(imprinted_endosperm_tem_genes$weighted_CG[which(imprinted_endosperm_tem_genes$class=="biparental")])

#mean(imprinted_endosperm_tem_genes$weighted_CHG[which(imprinted_endosperm_tem_genes$class=="maternal")])
#mean(imprinted_endosperm_tem_genes$weighted_CHG[which(imprinted_endosperm_tem_genes$class=="biparental")])

#endo teM genes prediction
endo_teM<-weighted_methylation%>%
  filter(gene%in%all_BM_5FC_core$V1)%>%
  mutate(class="endo teM genes")%>%
  select("gene","weighted_CG","weighted_CHG","class")


#fmegs
MEG_FMEGs <- read.delim("~/Documents/lab/megpeg/B73_MEGS_FMEGS.txt", header=T)

fmegs<-weighted_methylation%>%
  filter(gene%in%MEG_FMEGs$x)%>%
  mutate(class="FMEGs")%>%
  select("gene","weighted_CG","weighted_CHG","class")

#core_genes
core_gene<-read.csv("~/Documents/lab/b73_core_genes.txt", sep="")
core_genes_methylation<-weighted_methylation%>%
  filter(gene%in%core_gene$x)%>%
  mutate(class="Core Gene")%>%
  select("gene","weighted_CG","weighted_CHG","class")

all_methylation<-rbind(pollen_teM_genes_nm,imprinted_endosperm_tem_genes,endo_teM,fmegs,core_genes_methylation)

library(dplyr)
library(tidyr)
library(ggplot2)


df_long <- all_methylation %>%
  pivot_longer(cols = c(weighted_CG, weighted_CHG),
               names_to = "context", values_to = "methylation")%>%
  mutate(class=factor(class,levels=c("maternal","biparental","endo teM genes","DEG","non_DEG","FMEGs","Core Gene")))

library(ggplot2)
library(ggpubr)
ggplot(df_long, aes(x = methylation, y = class, fill = class)) +
  geom_boxplot(
    aes(color = class), # outline by class color
    fill = NA  ,
    outlier.shape = NA
  ) +
  # overlay points for classes other than "Core Gene"
  
  geom_jitter(
    data = subset(df_long, class != "Core Gene"),
    aes(color = class),
    width = 0, height = 0.3, size = 1.5, alpha = 0.6
  ) +
  facet_wrap(~context, scales = "free_y",
              labeller = as_labeller(c(
                "weighted_CG" = "CG ",
                "weighted_CHG" = "CHG "
              ))) +
  theme_bw(base_size = 16) +
  #scale_x_continuous(limits = c(0, 0.07)) +
  scale_x_continuous(limits = c(0, 0.12)) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey90")
  ) +
  labs(
    y = " ",
    x = "N-normalized methylation in TSS_centered 200bp",
    title = " "
  )+
  scale_y_discrete(limits = rev) +
  stat_compare_means(
    comparisons = list(c("maternal","biparental"),
                       c("DEG", "non_DEG")), 
    # which groups to compare
    method = "wilcox.test",
    label = "p.format", 
    method.args = list(alternative = "greater") ,
    #label.y = c(0.05, 0.05),
    label.y = c(0.09, 0.09),
    format.args = list(digits = 3,nsmall = 3)
     
  ) 
