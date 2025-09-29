#meta plot of CG, CHG, and CHH nucleotide frequency for BMEGs, FMEGs, and core genes
module load BEDTools/2.30.0-GCC-12.2.0

#using the gene flanking regions information in bed file, extract the sequences
bedtools getfasta -fi Zm-B73-REFERENCE-NAM-5.0.fa -bed tss_regions.bed -s -name > TSS.fa
bedtools getfasta -fi Zm-B73-REFERENCE-NAM-5.0.fa -bed tts_regions.bed -s -name > TTS.fa
