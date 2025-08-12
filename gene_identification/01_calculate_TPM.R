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

write.table(all_tpm,"~/Documents/lab/all_tpm.txt",quote=F)
