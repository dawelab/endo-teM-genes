library(dplyr)
b73_nodup_endo6_38_counts <- read.delim("~/Documents/lab/endo6_38_counts_dec.txt", comment.char="#")
Zm.B73.REFERENCE <- read.delim("~/Documents/lab/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3", header=FALSE, comment.char="#")

canonical.REFERENCE <- Zm.B73.REFERENCE %>% 
  filter(grepl("canonical_transcript=1", V9))

search_strings <-substring(canonical.REFERENCE$V9,4,23)

filtered_b73_endo6_38<-b73_nodup_endo6_38_counts[which(b73_nodup_endo6_38_counts$Geneid%in%search_strings),]
filtered_b73_endo6_38<-filtered_b73_endo6_38[,c(1:6,order(colnames(filtered_b73_endo6_38)[7:27])+6)]
rownames(filtered_b73_endo6_38)<-filtered_b73_endo6_38$Geneid
data<-filtered_b73_endo6_38[7:27]
meancounts<-filtered_b73_endo6_38[,1:6]
meancounts$D6<-data[, 1]
meancounts$D8<-data[, 2]
meancounts$D10<-data[, 3]
meancounts$D12<-rowMeans(data[, 4:5])
meancounts$D14<-rowMeans(data[, 6:7])
meancounts$D16<-data[, 8]
meancounts$D18<-rowMeans(data[, 9:10])
meancounts$D20<-data[, 11]
meancounts$D22<-data[, 12]
meancounts$D24<-rowMeans(data[, 13:14])
meancounts$D26<-data[, 15]
meancounts$D28<-data[, 16]
meancounts$D30<-data[, 17]
meancounts$D32<-data[, 18]
meancounts$D34<-data[, 19]
meancounts$D36<-data[, 20]
meancounts$D38<-data[, 21]
tpm <- function(x) (x/meancounts$Length)/(sum(x/meancounts$Length))*1e6
all_tpm <- lapply(meancounts[7:23], tpm)
all_tpm<-as.data.frame(all_tpm)
rownames(all_tpm) <- substr(rownames(meancounts), 1, nchar(rownames(meancounts)) - 5)
all_tpm$gene<-rownames(all_tpm)
B73.all <- read.csv("~/Documents/lab/megpeg/B73.all.csv")

b73.all.6_38<-merge(B73.all,all_tpm,by="gene",all=T)
b73.all.6_38 <- b73.all.6_38 %>%
  mutate(epiallele2 = case_when(
    mCG <= 0.05 & mCHG <= 0.05  ~  "UM",
    mCG >= 0.2 & mCHG <= 0.05  ~  "gbM",
    mCG >= 0.4 & mCHG >= 0.4 & cCHG >= 10 ~  "teM",
    TRUE ~ "ambiguous"         
  ))


#make figures by body-methylated groups############
all_BM_5FC_core <- read.table("~/Documents/lab/megpeg/all_BM_5FC_core.bed", quote="\"", comment.char="")
library(ggplot2)
library(tidyr)
library(dplyr)
bm_counts<-merge(all_BM_5FC_core,b73.all.6_38,by.x = "V1", by.y = "gene")
bm_counts[,c(25:41)]<-t(apply(bm_counts[,c(25:41)],1,function(x) (x - mean(x)) / sd(x)))
bm_counts<-bm_counts[,c(1,25:41)]

df <- bm_counts %>%
  pivot_longer(cols = -V1, names_to = "DAP", values_to = "TPM")%>%
  mutate(DAP = factor(DAP, levels = paste0("D", seq(6, 38, by = 2))))

ggplot(df, aes(x=DAP, y=TPM, group=V1)) +
  geom_line()+
  geom_point()

library(factoextra)

# Remove gene names for clustering
data_for_kmeans <- bm_counts[, 2:18]  # Only use normalized TPM values

# Find optimal k using the Elbow Method
fviz_nbclust(data_for_kmeans, kmeans, method = "wss") +
  ggtitle("Elbow Method for Optimal k")

fviz_nbclust(data_for_kmeans, kmeans, method = "silhouette") +
  ggtitle("Silhouette Method for Optimal k")

kmeans_result <- kmeans(bm_counts[,-1], centers = 2, nstart = 25)

bm_counts$Cluster <- as.factor(kmeans_result$cluster)
c1<-bm_counts[which(bm_counts$Cluster==1),]
c2<-bm_counts[which(bm_counts$Cluster==2),]

maxtpm<-b73.all.6_38
maxtpm$max<-apply(maxtpm[, c(25:41)],1, max)
c1<-merge(c1,maxtpm[,c(1,43)], by.x = "V1", by.y= "gene")
c2<-merge(c2,maxtpm[,c(1,43)], by.x = "V1", by.y= "gene")

#write.csv(c1$V1, "~/Documents/lab/megpeg/Re_ Mdr1 manuscript/cluster1.bed",quote = F,row.names = F)
#write.csv(c2$V1, "~/Documents/lab/megpeg/Re_ Mdr1 manuscript/cluster2.bed",quote = F,row.names = F)

df <- c1 %>%
  pivot_longer(cols = 2:18, names_to = "DAP", values_to = "TPM")%>%
  mutate(DAP = factor(DAP, levels = paste0("D", seq(6, 38, by = 2))))
df <- c2 %>%
  pivot_longer(cols = 2:18, names_to = "DAP", values_to = "TPM")%>%
  mutate(DAP = factor(DAP, levels = paste0("D", seq(6, 38, by = 2))))



ggplot(df, aes(x = DAP, y = TPM, group = V1, color = log2(max))) +  
  geom_line() +
  labs(x = "Days after pollination", y = "Z-score", color = "TPM") + 
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  stat_summary(aes(group = Cluster), fun = mean, geom = "line", size = 1.5, color = "black") +
  scale_color_gradient2(
    low = "blue", 
   mid = "grey",
    high = "red", 
    midpoint = log2(1000),
    #labels = function(x) round(2^(x))+
    breaks = log2(c( 10, 100, 1000,10000,100000)),
    labels = c( "10", "100", "1000","10000","100000"),
    limits = c(log2(10), log2(100000))
  )




df_long <- bm_counts %>%
  pivot_longer(cols = 2:18, names_to = "TimePoint", values_to = "TPM")%>%
  mutate(TimePoint = factor(TimePoint, levels = paste0("D", seq(6, 38, by = 2))))

ggplot(df_long, aes(x = TimePoint, y = TPM, group = V1, color = Cluster)) +
  geom_line(alpha = 0.4) +  # Faint lines for all genes
  stat_summary(aes(group = Cluster), fun = mean, geom = "line", size = 1.5, color = "black") +
  theme_minimal() +
  labs(x = "Days after pollination", y = "Z-score") +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


