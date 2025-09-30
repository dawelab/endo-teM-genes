#quantify number of exons (Figure 3C)
Zm.B73.REFERENCE <- read.delim("~/Documents/lab/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3", header=FALSE, comment.char="#")
mrna<-Zm.B73.REFERENCE[which(Zm.B73.REFERENCE$V3=="mRNA"),]
mrna<-mrna %>% 
  filter(grepl("canonical_transcript=1", V9))
mrna$gene<-substr(mrna$V9,4,23)
exons<-Zm.B73.REFERENCE[which(Zm.B73.REFERENCE$V3=="exon"),]
exons$gene<-substr(exons$V9,8,27)
exons <- exons %>%
  filter(gene%in%mrna$gene)
nexons<-as.data.frame(table(exons$gene))
nexons$genes<-substr(nexons$Var1,1,15)

core_gene <- read.csv("~/Documents/lab/b73_core_genes.txt", sep="")
endosperm_teM_genes <- read.table("~/Documents/lab/endosperm_teM_genes.txt", quote="\"", comment.char="")
FMEGS <- read.csv("~/Documents/lab/B73_MEGS_FMEGS.txt", sep="")
specific<- read.csv("~/Documents/lab/endosperm_specific_genes.txt", sep="")

coreexons<-nexons %>%
  filter(genes%in%core_gene$x)
coreexons$group<-rep("core_genes",nrow(coreexons))

specificexons<-nexons %>%
  filter(genes%in%specific$x)
specificexons$group<-rep("endosperm_specific_genes",nrow(specificexons))

nexons <- nexons %>%
  mutate(group = case_when(
    genes%in%endosperm_teM_genes$x  ~  "endosperm_teM_genes",
    genes%in%FMEGs$x  ~  "FMEGs",
    TRUE ~ "others"         
  ))

nexons<-nexons[-which(nexons$group=="others"),-1]
nexons<-rbind(nexons,coreexons[,-1],specificexons[,-1])

nexons$group_wrapped <- gsub("_", "\n", nexons$group)

ggplot(nexons, aes(x=group, y=Freq, fill=group)) +
  geom_boxplot(width=0.75)+
  theme_classic(base_size = 18)+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5)) +
  scale_y_continuous(trans = "log2",
                     breaks = c(1,2,5,10,20,40,80),
                     labels = c("1","2","5","10","20","40","80")) +
  labs(y="Number of exons",x="")+
  theme(legend.position = "none") +
  scale_x_discrete(
    limits = c("endosperm_teM_genes", "FMEGs", "endosperm_specific_genes", "core_genes"),
    labels = c("endosperm\nteM genes", "FMEGs", "endosperm\nspecific\ngenes", "core genes")
)



#length of elements (Figure 3D)

#generate intorn annotation file 
#ml  GenomeTools/1.6.2-GCC-11.3.0
#awk '{print "##sequence-region\t"$1"\t1\t"$2}' ./Zm-B73-REFERENCE-NAM-5.0.fa.fai > header.txt
#(head -n 1 Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3; cat header.txt; tail -n +2 Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3) > gff3_with_regions.gff3
#gt gff3 -addintrons -retainids -sortlines gff3_with_regions.gff3 > B73_introns.gff3



B73_introns <- read.table("~/Documents/lab/B73_introns.gff3", header=FALSE, quote="\"")

colnames(B73_introns)[9]<-"V9"
canonical.REFERENCE <- B73_introns %>% 
  filter(grepl("canonical_transcript=1", V9))

library(stringr)
canonical_ids <- str_extract(canonical.REFERENCE$V9, "Zm[0-9a-zA-Z_]+")

B73_introns$Parent <- str_extract(B73_introns$V9, "Parent=Zm[0-9a-zA-Z_]+") %>%
  str_replace("Parent=", "")
canonical_elements <- B73_introns %>%
  filter(Parent %in% canonical_ids)


canonical_elements$gene<-substring(canonical_elements$Parent,1,15)
canonical_elements$length<-canonical_elements$V5-canonical_elements$V4+1


BMEGS_length<-canonical_elements%>%
  filter(gene %in% endosperm_teM_genes$V1)%>%
  filter(V3 %in% c("intron","five_prime_UTR","three_prime_UTR","CDS")) %>%
  mutate(V3=case_when(V3 %in% c("five_prime_UTR", "three_prime_UTR")~"UTR",
                      TRUE ~ V3)) %>%
  group_by(gene, V3) %>%
  summarise(total_length = sum(length), .groups = "drop")

FMEGS_length<-canonical_elements%>%
  filter(gene %in% FMEGS$x)%>%
  filter(V3 %in% c("intron","five_prime_UTR","three_prime_UTR","CDS")) %>%
  mutate(V3=case_when(V3 %in% c("five_prime_UTR", "three_prime_UTR")~"UTR",
                      TRUE ~ V3)) %>%
  group_by(gene, V3) %>%
  summarise(total_length = sum(length), .groups = "drop")

core_gene_length<-canonical_elements%>%
  filter(gene %in% core_gene$x)%>%
  filter(V3 %in% c("intron","five_prime_UTR","three_prime_UTR","CDS")) %>%
  mutate(V3=case_when(V3 %in% c("five_prime_UTR", "three_prime_UTR")~"UTR",
                      TRUE ~ V3))%>%
  group_by(gene, V3) %>%
  summarise(total_length = sum(length), .groups = "drop")

BMEGS_length$geneset<-"Endo teM Genes"
FMEGS_length$geneset<-"FMEGs"
core_gene_length$geneset<-"Core Genes"

total_length<-bind_rows(BMEGS_length2,FMEGS_length2,core_gene_length2)
total_length$geneset <- factor(total_length$geneset, levels = c("Endo teM Genes", "FMEGs", "Core Genes"))
ggplot(total_length,aes(x=V3,y=(total_length/1000),fill = V3)) +
  geom_boxplot()+
  facet_wrap(~ geneset) +
  xlab("") +
  ylab("Element Length (kb)")+
  labs(fill = "element")+
  theme_bw(base_size = 16) +
  coord_cartesian(ylim = c(0, 5))+
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = 'top',
    legend.title = element_blank()
  )