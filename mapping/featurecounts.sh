module load Subread/2.0.6-GCC-11.3.0

featureCounts -p --countReadPairs  -t exon -g Parent -a /scratch/ys33815/sun/Zm-W22-REFERENCE-NRGENE-2.0_Zm00004b.1.gff3 \
-o ../nodup_endo6_38_counts.txt SRR1169626Aligned_nodup_q20.bam SRR1169633Aligned_nodup_q20.bam SRR1169641Aligned_nodup_q20.bam \
SRR1169627Aligned_nodup_q20.bam SRR1169634Aligned_nodup_q20.bam SRR1169642Aligned_nodup_q20.bam \
SRR1169628Aligned_nodup_q20.bam SRR1169635Aligned_nodup_q20.bam SRR1169643Aligned_nodup_q20.bam \
SRR1169629Aligned_nodup_q20.bam SRR1169636Aligned_nodup_q20.bam SRR1169644Aligned_nodup_q20.bam \
SRR1169630Aligned_nodup_q20.bam SRR1169637Aligned_nodup_q20.bam SRR1169645Aligned_nodup_q20.bam \
SRR1169631Aligned_nodup_q20.bam SRR1169639Aligned_nodup_q20.bam SRR1169646Aligned_nodup_q20.bam \
SRR1169632Aligned_nodup_q20.bam SRR1169640Aligned_nodup_q20.bam SRR1169647Aligned_nodup_q20.bam
