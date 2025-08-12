#genes are MEGs at  at least one time point when the expression is more than 1 rpm
#then remove the genes that are identified as both MEGs and PEGs
#transfer the B73 genes to W22
#select the genes that are overlapping with DMR (in gene body or 1kb flanking regions)
#transform W22 to B73
#select the core genes
#remove the teM endosperm genes

B73_W22_maternalpreference_with_avg_rpm <- read.csv("~/Documents/lab/B73_W22_maternalpreference_with_avg_rpm.csv")

B73_W22_megs<-B73_W22_maternalpreference_with_avg_rpm%>%
  mutate(exp11=case_when(pmax(A_mean11,B_mean11)<1 ~ NA_character_, TRUE ~ exp11),
         exp14=case_when(pmax(A_mean14,B_mean14)<1 ~ NA_character_, TRUE ~ exp14),
         exp17=case_when(pmax(A_mean17,B_mean17)<1 ~ NA_character_, TRUE ~ exp17),
         exp21=case_when(pmax(A_mean21,B_mean21)<1 ~ NA_character_, TRUE ~ exp21)
         )%>%
  filter(
    (exp11=="maternal"|exp14=="maternal"|exp17=="maternal"|exp21=="maternal") &
   !(exp11 == "paternal" | exp14 == "paternal" | exp17 == "paternal" | exp21 == "paternal"))

W22_B73_combine <- read.delim("~/Documents/lab/W22_B73_combine.csv", header=FALSE)

converted <- B73_W22_megs %>%
  mutate(genome=case_when(
    substring(feature,8,9)=="eb" ~ "B73",
    TRUE~"W22"
  ))

converted<-merge(converted,W22_B73_combine,by.x = "feature",by.y = "V11",all.x = T)
converted<-converted %>%
  mutate(W22=case_when(
    genome=="B73"~V10,
    genome=="W22"~feature
  ))

Zm.W22.REFERENCE.NRGENE.2.0_Zm00004b.1 <- read.delim("~/Documents/lab/Zm-W22-REFERENCE-NRGENE-2.0_Zm00004b.1.gff3", header=FALSE)  

converted_megs_regions<-Zm.W22.REFERENCE.NRGENE.2.0_Zm00004b.1%>%
  filter(V3=="gene"&substring(V9,4,17)%in%converted$W22 )


converted_megs_regions <- converted_megs_regions %>%
  select(V1, V4, V5, V6, V7, V8, V9) %>%                      # Keep selected columns
  mutate(
    V4 = pmax(V4 - 1000, 0),                                 # Subtract 1000, min 0
    V5 = V5 + 1000,                                          # Add 1000
    V1 = substring(V1, 4)                                    # Remove 'chr' prefix
  )

#write.table(converted_megs_regions,"~/Documents/lab/megpeg/meg_regions.bed",quote = F,row.names = F,
#                     col.names = F,sep = "\t")

#module load BEDTools/2.30.0-GCC-12.2.0
#bedtools intersect -a meg_regions.bed -b /scratch/ys33815/sun/wt_endo-embryo_hyperDMRs.bed -u > meg_DMR.bed

meg_DMR <- read.delim("~/Documents/lab/megpeg/meg_DMR.bed", header=FALSE)
meg_DMR$feature<-substring(meg_DMR$V7,4,17)

B73_megs_DMR<-merge(meg_DMR,W22_B73_combine,by.x="feature",by.y="V10")
B73.all <- read.csv("~/Documents/lab/megpeg/B73.all.csv")  

B73_megs_core<-merge(B73_megs_DMR,B73.all,by.x="V11",by.y="gene",all.x=T) %>%
  filter(class=="Core Gene")%>%
  select(c(1,2,19:43))

B73_megs_core <- B73_megs_core %>%
  mutate(epiallele2 = case_when(
    mCG <= 0.05 & mCHG <= 0.05  ~  "UM",
    mCG >= 0.2 & mCHG <= 0.05  ~  "gbM",
    mCG >= 0.4 & mCHG >= 0.4 & cCHG >= 10 ~  "teM",
    TRUE ~ "ambiguous"         
  ))


B73_megs_FMEGs<-B73_megs_core%>%
  filter(epiallele2!="teM")
  
  
  
