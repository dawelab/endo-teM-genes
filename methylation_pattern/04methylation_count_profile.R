library(Biostrings)
count_methylation_contexts <- function(name, sequence, bin_size = 100) {
  sequence <- toupper(sequence)  # Ensure uppercase
  n <- nchar(sequence)
  bins <- seq(1, n, by = bin_size)
  
  data.frame(
    gene = name,
    bin = seq_along(bins),
    start = bins,
    end = pmin(bins + bin_size - 1, n),
    segment = sapply(seq_along(bins), function(i) substr(sequence, bins[i], min(bins[i] + bin_size - 1, n)))
  ) %>%
    rowwise() %>%
    mutate(
      CG_percent  = str_count(segment, "CG"),
      CHG_percent = str_count(segment, "C[ATC]G"),
      CHH_percent = str_count(segment, "C[ATC][ATC]"),
    ) %>%
    ungroup() %>%
    select(gene, bin, CG_percent, CHG_percent, CHH_percent)
}

count_methylation_contexts2 <- function(name, sequence, bin_size = 100) {
  sequence <- toupper(sequence)
  n <- nchar(sequence)
  bins <- seq(1, n, by = bin_size)
  
  bin_data <- data.frame(
    gene = name,
    bin = seq_along(bins),
    start = bins,
    end = pmin(bins + bin_size - 1, n),
    segment = sapply(seq_along(bins), function(i) substr(sequence, bins[i], min(bins[i] + bin_size - 1, n)))
  )
  
  bin_data <- bin_data %>%
    rowwise() %>%
    mutate(
      CG_percent = {
        s <- segment
        sum(sapply(1:(nchar(s) - 1), function(i) substr(s, i, i + 1) == "CG"))
      },
      CHG_percent = {
        s <- segment
        sum(sapply(1:(nchar(s) - 2), function(i) {
          b1 <- substr(s, i, i)
          b2 <- substr(s, i + 1, i + 1)
          b3 <- substr(s, i + 2, i + 2)
          b1 == "C" && b2 %in% c("A", "T", "C") && b3 == "G"
        }))
      },
      CHH_percent = {
        s <- segment
        sum(sapply(1:(nchar(s) - 2), function(i) {
          b1 <- substr(s, i, i)
          b2 <- substr(s, i + 1, i + 1)
          b3 <- substr(s, i + 2, i + 2)
          b1 == "C" && b2 %in% c("A", "T", "C") && b3 %in% c("A", "T", "C")
        }))
      }
    ) %>%
    ungroup() %>%
    select(gene, bin, CG_percent, CHG_percent, CHH_percent)
  
  return(bin_data)
}
tss <- readDNAStringSet("~/Documents/lab/megpeg/TSS2.fa", format = "fasta")
seq_names <- names(tss)
seqs <- as.character(tss)
all_results_tss <- purrr::map2_dfr(seq_names, seqs, count_methylation_contexts2)

tts <- readDNAStringSet("~/Documents/lab/megpeg/TTS2.fa", format = "fasta")
seq_names <- names(tts)
seqs <- as.character(tts)
all_results_tts <- purrr::map2_dfr(seq_names, seqs, count_methylation_contexts2)

all_results_tss$region <- "TSS"
all_results_tts$region <- "TTS"

tss_profile <- all_results_tss %>%
  mutate(position = (bin - 31) * 100)  # from -1000 to +900

tts_profile <- all_results_tts %>%
  mutate(position = (bin - 1) * 100 + 2000)  # from +2000 to +3900 (adds a visual gap)

combined_profile <- bind_rows(tss_profile, tts_profile)

matched_rows <- vector("list", 40*length(all_BM_5FC_core$V1))

for (i in seq_along(all_BM_5FC_core$V1)) {
  gene <- all_BM_5FC_core$V1[i]
  matches <- grepl(gene,combined_profile$gene, fixed = TRUE)
  if (any(matches)) {
    matched_rows[[i]] <- combined_profile[matches, ]
  }
}

bmegs_tss_tts<- do.call(rbind, matched_rows)
bmegs_tss_tts$sample<-"BMEG"

matched_rows <- vector("list", 40*length(MEG_FMEG$V1))

for (i in seq_along(MEG_FMEG$V1)) {
  gene <- MEG_FMEG$V1[i]
  matches <- grepl(gene,combined_profile$gene, fixed = TRUE)
  if (any(matches)) {
    matched_rows[[i]] <- combined_profile[matches, ]
  }
}

fmegs_tss_tts<- do.call(rbind, matched_rows)
fmegs_tss_tts$sample<-"FMEG"

matched_rows <- vector("list", 40*length(core_gene$V1))

for (i in seq_along(core_gene$V1)) {
  gene <- core_gene$V1[i]
  matches <- grepl(gene,combined_profile$gene, fixed = TRUE)
  if (any(matches)) {
    matched_rows[[i]] <- combined_profile[matches, ]
  }
}

coregene_tss_tts<- do.call(rbind, matched_rows)
coregene_tss_tts$sample<-"core gene"

matched_rows<-vector("list", 40*length(mpg$`gene ID`))
for (i in seq_along(mpg$`gene ID`)) {
  gene<-mpg$`gene ID`[i]
  matches<-grepl(gene, combined_profile$gene, fixed = T)
  if (any(matches)) {
    matched_rows[[i]]<-combined_profile[matches,]
  }
}

mpg_tss_tts<-do.call(rbind,matched_rows)
mpg_tss_tts$sample<-"MPG"

combined<-rbind(bmegs_tss_tts,fmegs_tss_tts,mpg_tss_tts,coregene_tss_tts)  

meta_profile <- combined %>%
  pivot_longer(cols = c("CG_percent", "CHG_percent", "CHH_percent"),
               names_to = "Context", values_to = "Percent") %>%
  group_by(position, Context, region, sample) %>%
  summarise(Mean = mean(Percent, na.rm = TRUE), .groups = "drop") %>%
  mutate(group = paste(sample, region, Context, sep = "_"))

meta_profile$sample <- factor(meta_profile$sample, levels = c("BMEG","FMEG","MPG","core gene"))
#Total bases: 2182075994
#CG: 93930955(4.3%)
#CHG: 81053743(3.71%)
#CHH: 325052331(14.89%)
ggplot(filter(meta_profile, Context == "CHH_percent"), aes(x = position, y = Mean, color = sample, group = group)) +
  geom_line(size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 3000, linetype = "dashed", color = "gray") +
  geom_hline(yintercept = 14.89, linetype = "dashed", color = "black")+
  scale_x_continuous(
    breaks = c(-3000, 0, 3000, 6000),
    labels = c("-3 kb", "TSS", "TTS", "+3 kb")
  ) +
  labs(
    x = "Position relative to TSS–TTS (bp)",
    y = "Mean % cytosine context",
    title = " CHH around TSS and TTS"
  ) +
  theme_minimal() +
  theme(text = element_text(size = 20))
                     
           
