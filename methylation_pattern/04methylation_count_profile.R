library(Biostrings)
library(dplyr)

count_methylation_contexts <- function(name, sequence, bin_size = 100) {
  sequence <- toupper(sequence)
  n <- nchar(sequence)
  bins <- seq(1, n, by = bin_size)
  
  bin_data <- data.frame(
    gene = name,
    bin = seq_along(bins),
    start = bins,
    end = pmin(bins + bin_size -1, n),
    segment = sapply(seq_along(bins), function(i) substr(sequence, bins[i], min(bins[i] + bin_size - 1, n)))
  )
  
  bin_data <- bin_data %>%
    rowwise() %>%
    mutate(
      CG_percent = {
        s1 <- segment
        s2 <- as.character(reverseComplement(DNAString(s1)))
        sum(sapply(c(s1, s2), function(s) {
          sum(sapply(1:(nchar(s) - 1), function(i) substr(s, i, i + 1) == "CG"))
        }))
      },
      CHG_percent = {
        s1 <- segment
        s2 <- as.character(reverseComplement(DNAString(s1)))
        sum(sapply(c(s1, s2), function(s) {
          sum(sapply(1:(nchar(s) - 2), function(i) {
            b1 <- substr(s, i, i)
            b2 <- substr(s, i + 1, i + 1)
            b3 <- substr(s, i + 2, i + 2)
            b1 == "C" && b2 %in% c("A", "T", "C") && b3 == "G"
          }))
        }))
      },
      CHH_percent = {
        s1 <- segment
        s2 <- as.character(reverseComplement(DNAString(s1)))
        sum(sapply(c(s1, s2), function(s) {
          sum(sapply(1:(nchar(s) - 2), function(i){
            b1 <- substr(s, i, i)
            b2 <- substr(s, i + 1, i + 1)
            b3 <- substr(s, i + 2, i + 2)
            b1 == "C" && b2 %in% c("A", "T", "C") && b3 %in% c("A", "T", "C")
          }))
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
all_results_tss <- purrr::map2_dfr(seq_names, seqs, count_methylation_contexts)

tts <- readDNAStringSet("~/Documents/lab/megpeg/TTS2.fa", format = "fasta")
seq_names <- names(tts)
seqs <- as.character(tts)
all_results_tts <- purrr::map2_dfr(seq_names, seqs, count_methylation_contexts)

all_results_tss$region <- "TSS"
all_results_tts$region <- "TTS"

tss_profile <- all_results_tss %>%
  mutate(position = (bin - 31) * 100)  # from -3000 to +900

tts_profile <- all_results_tts %>%
  mutate(position = (bin - 1) * 100 + 2000)  # from +2000 to +5900

combined_profile <- bind_rows(tss_profile,tts_profile)

write.table(combined_profile,"~/Documents/lab/megpeg/tss_tts_counts_doublestrand.csv",quote = F,
            col.names = T,row.names = F,sep = "\t")
                     
           
