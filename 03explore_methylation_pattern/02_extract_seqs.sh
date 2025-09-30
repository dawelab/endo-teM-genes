module load BEDTools/2.30.0-GCC-12.2.0

#using the gene flanking regions information in bed file, extract the sequences
bedtools getfasta -fi Zm-B73-REFERENCE-NAM-5.0.fa -bed tss_regions_200bp.bed -s -name > TSS_200bp.fa

#extract the sequences in the gene body regions
bedtools getfasta -fi /scratch/ys33815/Zm-B73-REFERENCE-NAM-5.0.fa -bed gene_body_regions_ALL.bed -s -name > gene_body.fa
