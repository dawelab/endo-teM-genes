#check the tissue specificity(Figure 3A)

#FMEGs
FMEGS <- read.delim("~/Documents/lab/megpeg/B73_MEGS_FMEGS.txt", header=FALSE)
FMEGS <-FMEGS[,1]
FMEGS_FC<-b73.all.6_38[which(b73.all.6_38$gene%in%FMEGS),]
FMEGS_FC$maxendo<-apply(FMEGS_FC[,25:41],1,max)#find the max endosperm expression
FMEGS_FC$max<-apply(FMEGS_FC[,c(14,16:23)],1,max)#find the max other tissues expression
FMEGS_FC$FC<-FMEGS_FC$maxendo/(FMEGS_FC$max+1)

#endosperm teM genes
endo_tem <- read.table("~/Documents/lab/endosperm_teM_genes.txt", quote="\"", comment.char="")
endo_tem_FC<-b73.all.6_38[which(b73.all.6_38$gene%in%endo_tem$V1),]
endo_tem_FC$maxendo<-apply(endo_tem_FC[,25:41],1,max)
endo_tem_FC$max<-apply(endo_tem_FC[,c(14,16:23)],1,max)
endo_tem_FC$FC<-(endo_tem_FC$maxendo+1)/(endo_tem_FC$max+1)

#endosperm specific genes
specific<- read.csv("~/Documents/lab/endosperm_specific_genes.txt", sep="")
specific<-merge(specific,b73.all.6_38, by.x="x",by.y="gene")
specific$maxendo<-apply(specific[,25:41],1,max)
specific$max<-apply(specific[,c(14,16:23)],1,max)
specific$FC<-(specific$maxendo+1)/(specific$max+1)
specific<-specific[which(specific$maxendo>=20&specific$FC>=5),]

#core genes
core_gene<-read.csv("~/Documents/lab/b73_core_genes.txt", sep="")
core_gene<-merge(core_gene,b73.all.6_38,by.x="x",by.y="gene")
core_gene$maxendo<-apply(core_gene[,25:41],1,max)
core_gene$max<-apply(core_gene[,c(14,16:23)],1,max)
core_gene$FC<-(core_gene$maxendo+1)/(core_gene$max+1)
  
fc1<-as.data.frame(endo_tem_FC$FC)
colnames(fc1)<-"FC"
fc1$group<-rep("endosperm_teM_genes",nrow(fc1))
fc2<-as.data.frame(FMEGS_FC$FC)
colnames(fc2)<-"FC"
fc2$group<-rep("FMEGs",nrow(fc2))
fc3<-as.data.frame(specific$FC)
colnames(fc3)<-"FC"
fc3$group<-rep("endosperm_specific_genes",nrow(fc3))
fc4<-as.data.frame(core_gene$FC)
colnames(fc4)<-"FC"
fc4$group<-rep("core_genes",nrow(fc4))
FC<-rbind(fc1,fc2,fc3,fc4)

FC$group_wrapped <- gsub("_", "\n", FC$group)


ggplot(FC, aes(x=group, y=FC+1, fill=group)) +
  geom_boxplot(width = 0.75)+
  scale_y_continuous(trans = "log2",
                     breaks =c(1,10,100,1000,10000,100000),
                     labels = c("1","10","100","1000","10000","100000"))+
  theme_classic(base_size = 18,)+
  theme(
    legend.position = "none") +
  labs(y="Fold Change \n endosperm/sporophyte",x="")+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5)) +
  scale_x_discrete(
    limits = c("endosperm_teM_genes", "FMEGs", "endosperm_specific_genes", "core_genes"),
    labels = c("endosperm\nteM genes", "FMEGs", "endosperm\nspecific\ngenes", "core genes")
)



#check the max expression level of each gene set(Figure 3B)

fc1<-as.data.frame(en_tem_FC$maxendo)
colnames(fc1)<-"TPM"
fc1$group<-rep("endosperm_teM_genes",nrow(fc1))
fc2<-as.data.frame(FMEGS_FC$maxendo)
colnames(fc2)<-"TPM"
fc2$group<-rep("FMEGs",nrow(fc2))

fc3<-as.data.frame(specific$maxendo)
colnames(fc3)<-"TPM"
fc3$group<-rep("endosperm_specific_genes",nrow(fc3))
fc4<-as.data.frame(core_gene$maxendo)
colnames(fc4)<-"TPM"
fc4$group<-rep("core_genes",nrow(fc4))
TPM<-rbind(fc1,fc2,fc3,fc4)

TPM$group_wrapped <- gsub("_", "\n", FC$group)


ggplot(TPM, aes(x=group, y=TPM+1, fill=group)) +
  geom_boxplot(width = 0.75)+
  scale_y_continuous(trans = "log2",
                     breaks =c(1,10,100,1000,10000,100000),
                     labels = c("1","10","100","1000","10000","100000"))+
  theme_classic(base_size = 18,)+
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5)) +
  labs(y="TPM in endosperm",x="")+
  theme(legend.position = "none") +
  scale_x_discrete(
    limits = c("endosperm_teM_genes", "FMEGs", "endosperm_specific_genes", "core_genes"),
    labels = c("endosperm\nteM genes", "FMEGs", "endosperm\nspecific\ngenes", "core genes")
)


#stacked bar chart of endosperm tem genes expreesion(Figure 3E)
b73_core_genes <- read.csv("~/Documents/lab/b73_core_genes.txt", sep="")
endosperm_teM_genes <- read.table("~/Documents/lab/endosperm_teM_genes.txt", quote="\"", comment.char="")
FMEGS <- read.csv("~/Documents/lab/B73_MEGS_FMEGS.txt", sep="")

stacked_data<-b73.all.6_38%>%
  filter(gene%in%b73_core_genes$x)%>%
  mutate(group=ifelse(gene%in%endosperm_teM_genes$V1,"Endo teM Genes","Other Core Genes"))%>%
  mutate(group=ifelse(gene%in%FMEGS$x,"FMEGs",group))

df_long<-stacked_data%>%
  pivot_longer(cols = 25:41, names_to = "TimePoint", values_to = "TPM")%>%
  mutate(TimePoint = factor(substring(TimePoint,2), levels = seq(6, 38, by = 2)))


ggplot(df_long, aes(x = TimePoint, y = TPM, fill = group)) +
  geom_bar(stat = "identity",position = "fill") +
  #scale_y_continuous(labels = scales::percent_format()) + 
  scale_y_continuous(labels = scales::percent_format(scale = 100, suffix = "", accuracy = 1))+
  scale_fill_manual(values = c("Endo teM Genes" = "orange","FMEGs"="black", "Other Core Genes" = "lightblue")) +
  theme_minimal() +
  labs(title = "TPM composition at different endosperm stages",
       x = "Days After Pollination (DAP)", y = "Total TPM") +
  theme(text = element_text(size = 18))+
  theme(#legend.position="top",
                #legend.box = "horizontal",
                 legend.title = element_blank()
)
