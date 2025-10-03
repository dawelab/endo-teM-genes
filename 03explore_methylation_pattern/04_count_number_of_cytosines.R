library(Biostrings)
library(dplyr)

count_methylation <- function(name, sequence) {
  sequence <- toupper(sequence)
  seq_rc <- as.character(reverseComplement(DNAString(sequence)))
  gene_start <- 3
  gene_end   <- nchar(sequence) - 2
  count_CG <- function(seq) {
    sum(sapply(gene_start:gene_end, function(i) {
      if (i <= nchar(seq) - 1) {
        substr(seq, i, i) == "C" && substr(seq, i + 1, i + 1) == "G"
      } else FALSE
    }))
  }
  
  count_CHG <- function(seq) {
    sum(sapply(gene_start:gene_end, function(i) {
      if (i <= nchar(seq) - 2) {
        b1 <- substr(seq, i, i)
        b2 <- substr(seq, i + 1, i + 1)
        b3 <- substr(seq, i + 2, i + 2)
        b1 == "C" && b2 %in% c("A","T","C") && b3 == "G"
      } else FALSE
    }))
  }
  
  count_CHH <- function(seq) {
    sum(sapply(gene_start:gene_end, function(i) {
      if (i <= nchar(seq) - 2) {
        b1 <- substr(seq, i, i)
        b2 <- substr(seq, i + 1, i + 1)
        b3 <- substr(seq, i + 2, i + 2)
        b1 == "C" && b2 %in% c("A","T","C") && b3 %in% c("A","T","C")
      } else FALSE
    }))
  }
  
  
  CG_total  <- count_CG(sequence)  + count_CG(seq_rc)
  CHG_total <- count_CHG(sequence) + count_CHG(seq_rc)
  CHH_total <- count_CHH(sequence) + count_CHH(seq_rc)
  
  tibble(
    gene = name,
    length_bp = nchar(sequence)-4,
    CG_count  = CG_total,
    CHG_count = CHG_total,
    CHH_count = CHH_total
  )
}


#count cytosines in 200bp regions centered TSS
TSS_200bp <- readDNAStringSet("./TSS_200bp.fa", format = "fasta")
seq_names <- names(TSS_200bp)
seqs <- as.character(TSS_200bp)
all_results_TSS_200bp<- purrr::map2_dfr(seq_names, seqs, count_methylation )


write.table(all_results_TSS_200bp"./tss_counts_doublestrand_200bp.csv",quote = F,
            col.names = T,row.names = F,sep = "\t")


           
#count cytosines in gene body regions
gene_body <- readDNAStringSet("./gene_body.fa", format = "fasta")
seq_names <- names(gene_body)
seqs <- as.character(gene_body)
all_results_gene_body<- purrr::map2_dfr(seq_names, seqs, count_methylation )


write.table(all_results_gene_body,"./counts_doublestrand_genebody.csv",quote = F,
            col.names = T,row.names = F,sep = "\t")
