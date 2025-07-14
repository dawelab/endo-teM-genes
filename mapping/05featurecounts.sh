module load Subread/2.0.6-GCC-11.3.0


featureCounts -p --countReadPairs -t exon -O -g Parent -a /scratch/ys33815/endo6_38/b73bam/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 \
-o ../endo6_38_counts_dec.txt SRR1169626Aligned.sortedByCoord.out.bam SRR1169633Aligned.sortedByCoord.out.bam SRR1169641Aligned.sortedByCoord.out.bam \
SRR1169627Aligned.sortedByCoord.out.bam SRR1169634Aligned.sortedByCoord.out.bam SRR1169642Aligned.sortedByCoord.out.bam \
SRR1169628Aligned.sortedByCoord.out.bam SRR1169635Aligned.sortedByCoord.out.bam SRR1169643Aligned.sortedByCoord.out.bam \
SRR1169629Aligned.sortedByCoord.out.bam SRR1169636Aligned.sortedByCoord.out.bam SRR1169644Aligned.sortedByCoord.out.bam \
SRR1169630Aligned.sortedByCoord.out.bam SRR1169637Aligned.sortedByCoord.out.bam SRR1169645Aligned.sortedByCoord.out.bam \
SRR1169631Aligned.sortedByCoord.out.bam SRR1169639Aligned.sortedByCoord.out.bam SRR1169646Aligned.sortedByCoord.out.bam \
SRR1169632Aligned.sortedByCoord.out.bam SRR1169640Aligned.sortedByCoord.out.bam SRR1169647Aligned.sortedByCoord.out.bam 
