# substract the 200bp regions centered TSS 
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
    TSS_start = ifelse(V7 == "+", V4 - 100-1-2, V5 - 100-1-2),#bed file is 0-based
    TSS_end   = ifelse(V7 == "+", V4 + 100+2, V5 + 100+2),# need 2bp extension to deal with edge cytosines
    TSS_start = ifelse(TSS_start < 0, 0, TSS_start)  
  )

TTS <- canonical.mrna %>%
  mutate(
    TTS_start = ifelse(V7 == "+", V5 - 100-1-2, V4 - 100-1-2),#bed file is 0-based
    TTS_end   = ifelse(V7 == "+", V5 + 100+2, V4 + 100+2), # need 2bp extension to deal with edge cytosines
    TTS_start = ifelse(TTS_start < 0, 0, TTS_start)  
  )

TSS_regions<-TSS[,c(1,11,12,10,6:8)]


#TTS_regions<-TTS[,c(1,11,12,10,6:8)]


write.table(TSS_regions, "~/Documents/lab/megpeg/tss_regions_200bp.bed",quote = F,row.names = F,
            col.names = F,sep = "\t")
#write.table(TTS_regions, "~/Documents/lab/megpeg/tts_regions_200bp.bed",quote = F,row.names = F,
#            col.names = F,sep = "\t")

