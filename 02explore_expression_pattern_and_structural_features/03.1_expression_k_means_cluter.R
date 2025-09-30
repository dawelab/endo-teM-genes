b73.all.6_38 <- read.csv("~/Documents/lab/b73.all.6_38.txt", sep="")

b73.all.6_38 <- b73.all.6_38 %>%
  mutate(epiallele2 = case_when(
    mCG <= 0.05 & mCHG <= 0.05  ~  "UM",
    mCG >= 0.2 & mCHG <= 0.05  ~  "gbM",
    mCG >= 0.4 & mCHG >= 0.4 & cCHG >= 10 ~  "teM",
    TRUE ~ "ambiguous"         
  ))
  

#make figures by gene clusters############
all_BM_5FC_core <- read.table("~/Documents/lab/endosperm_teM_genes.txt", quote="\"", comment.char="")
library(ggplot2)
library(tidyr)
library(dplyr)
bm_counts<-merge(all_BM_5FC_core,b73.all.6_38,by.x = "V1", by.y = "gene")
#calculate Z-score
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

# Find optimal k for K-means cluster
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

df1 <- c1 %>%
  pivot_longer(cols = 2:18, names_to = "DAP", values_to = "TPM")%>%
  mutate(DAP = factor(DAP, levels = paste0("D", seq(6, 38, by = 2))))
df2 <- c2 %>%
  pivot_longer(cols = 2:18, names_to = "DAP", values_to = "TPM")%>%
  mutate(DAP = factor(DAP, levels = paste0("D", seq(6, 38, by = 2))))



ggplot(df1, aes(x = DAP, y = TPM, group = V1, color = log2(max))) +  
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



#Select several genes in each cluster showing the real TPM 
  set<-c("Zm00001eb302040","Zm00001eb005200",
       "Zm00001eb375260","Zm00001eb042960",
       "Zm00001eb094330", "Zm00001eb166950",
       "Zm00001eb170070", "Zm00001eb313790")

setcounts<-b73.all.6_38[which(b73.all.6_38$gene%in%set),c(1,25:41)]
setcounts<-setcounts %>%
  mutate(genename = case_when(
           gene=="Zm00001eb302040" ~ "meg1",
           gene=="Zm00001eb375260" ~ "myb",
           gene=="Zm00001eb042960" ~ "tcrr2",
           gene=="Zm00001eb094330" ~ "pbf",
           gene=="Zm00001eb166950" ~ "19kD α-zein",
           gene=="Zm00001eb170070" ~ "22 kDa α-zein",
           gene=="Zm00001eb313790" ~ "50 kDa γ-zein",
           gene=="Zm00001eb005200" ~ "esr2",
           TRUE ~ "none"
         ))

df <- setcounts %>%
  pivot_longer(cols = 2:18, names_to = "DAP", values_to = "TPM")%>%
  mutate(DAP = factor(substring(DAP,2,), levels = seq(6, 38, by = 2)))

df$genename <- factor(df$genename, levels = c("19 kDa α-zein", "22 kDa α-zein", 
                                        "50 kDa γ-zein", "pbf", 
                                        "meg1", "myb","tcrr2", "esr2" ))

df$Cluster <- ifelse(df$genename %in% c("19 kDa α-zein", "22 kDa α-zein", "50 kDa γ-zein", "pbf"), "Cluster 2", "Cluster 1")


ggplot(df, aes(x=DAP, y=TPM + 1, group=genename,color=genename)) +
  geom_line() +
  scale_y_continuous(trans = "log2", 
                     labels = c("1","10", "100", "1000","10000","100000"),
                     breaks = c(1,10, 100, 1000, 10000, 100000)) +
  theme_classic(base_size = 18) +
  labs(x="Days After Pollination (DAP)",y= "TPM")+
  theme(#legend.position="top",
        #legend.box = "horizontal",
        legend.title = element_blank())




