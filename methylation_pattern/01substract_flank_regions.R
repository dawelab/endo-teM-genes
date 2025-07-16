library(dplyr)
library(stringr)
Zm.B73.REFERENCE <- read.delim("~/Documents/lab/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3", header=FALSE, comment.char="#")
canonical.mrna <- Zm.B73.REFERENCE %>% 
  filter(grepl("canonical_transcript=1", V9) & V3=="mRNA")
canonical.mrna <- canonical.mrna %>%
  mutate(ID = str_extract(V9, "ID=[^;]+")) %>%
  mutate(ID = str_replace(ID, "ID=", ""))

TSS <- canonical.mrna %>%
  mutate(
    TSS_start = ifelse(V7 == "+", V4 - 3000, V5 - 1000),
    TSS_end   = ifelse(V7 == "+", V4 + 1000, V5 + 3000),
    TSS_start = ifelse(TSS_start < 0, 0, TSS_start)  
  )

TTS <- canonical.mrna %>%
  mutate(
    TTS_start = ifelse(V7 == "+", V5 - 1000, V4 - 3000),
    TTS_end   = ifelse(V7 == "+", V5 + 3000, V4 + 1000),
    TTS_start = ifelse(TTS_start < 0, 0, TTS_start)  
  )

TSS_regions<-TSS[,c(1,11,12,10,6:8)]


TTS_regions<-TTS[,c(1,11,12,10,6:8)]

#region.bed is 1 kb centered tss and tts, region2 is 3kb upstream
write.table(TSS_regions, "~/Documents/lab/megpeg/tss_regions2.bed",quote = F,row.names = F,
                                col.names = F,sep = "\t")
write.table(TTS_regions, "~/Documents/lab/megpeg/tts_regions2.bed",quote = F,row.names = F,
            col.names = F,sep = "\t")
