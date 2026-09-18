module load SAMtools/1.16.1-GCC-11.3.0

samtools merge -o m_endo_H3K56ac.bam m1-H3K56ac.q20.sorted.bam m2-H3K56ac.q20.sorted.bam m4-H3K56ac.q20.sorted.bam
samtools merge -o w_endo_H3K56ac.bam w1-H3K56ac.q20.sorted.bam w2-H3K56ac.q20.sorted.bam w3-H3K56ac.q20.sorted.bam
samtools merge -o w_leaf_H3K56ac.bam J625xJ625-H3K56ac.q20.sorted.bam J625xJ625-H3K56ac-rep2.q20.sorted.bam J643self-H3K56ac.q20.sorted.bam
samtools merge -o m_endo_H3K27me3.bam m1-H3K27me3.q20.sorted.bam m2-H3K27me3.q20.sorted.bam
samtools merge -o w_endo_H3K27me3.bam w1-H3K27me3.q20.sorted.bam w3-H3K27me3.q20.sorted.bam
samtools merge -o w_leaf_H3K27me3.bam J625xJ625-H3K27me3.q20.sorted.bam J625xJ625-H3K27me3-39157.q20.sorted.bam

ml deepTools/3.5.5-gfbf-2023a
ml SAMtools/1.21-GCC-13.3.0

for i in w_endo_H3K27me3 w_endo_H3K56ac w_leaf_H3K27me3 w_leaf_H3K56ac W22-6; do
    bamCoverage -b ../${i}.bam -o ${i}.bw -bs 10 --normalizeUsing RPGC --effectiveGenomeSize 2300000000 -p 10
done

bigwigCompare \
  -b1 w_endo_H3K27me3.bw \
  -b2 W22-6.bw \
  --operation ratio \
  --binSize 10 \
  -o w_endo_H3K27me3_ratio_w22.bw \
  -p 10


bigwigCompare \
  -b1 w_endo_H3K56ac.bw \
  -b2 W22-6.bw \
  --operation ratio \
  --binSize 10 \
  -o w_endo_H3K56ac_ratio_w22.bw \
  -p 10

bigwigCompare \
  -b1 w_leaf_H3K27me3.bw \
  -b2 W22-6.bw \
  --operation ratio \
  --binSize 10 \
  -o w_leaf_H3K27me3_ratio_w22.bw \
  -p 10

bigwigCompare \
  -b1 w_leaf_H3K56ac.bw \
  -b2 W22-6.bw \
  --operation ratio \
  --binSize 10 \
  -o w_leaf_H3K56ac_ratio_w22.bw \
  -p 10

bigwigCompare \
  -b1 w_leaf_H3K56ac.bw \
  -b2 W22-6.bw \
  --operation ratio \
  --binSize 10 \
  -o w_leaf_H3K56ac_ratio_w22.bw \
  -p 10

bigwigCompare \
  -b1 w2-IgG.q20.sorted.bw \
  -b2 W22-6.bw \
  --operation ratio \
  --binSize 10 \
  -o endo_IgG_ratio_w22.bw \
  -p 10

bigwigCompare \
  -b1 J625xJ625-IgG.q20.sorted.bw \
  -b2 W22-6.bw \
  --operation ratio \
  --binSize 10 \
  -o leaf_IgG_ratio_w22.bw \
  -p 10


computeMatrix reference-point \
  -S endo_IgG_ratio_w22.bw leaf_IgG_ratio_w22.bw w_endo_H3K27me3_ratio_w22.bw w_endo_H3K56ac_ratio_w22.bw w_leaf_H3K27me3_ratio_w22.bw w_leaf_H3K56ac_ratio_w22.bw \
  -R cluster1_regions.bed cluster2_regions.bed endosperm_teM_genes_regions.bed endosperm_specific_genes_regions.bed pollen_teM_genes_regions.bed core_genes_regions.bed FMEGs_region.bed\
  -p 10 \
  --referencePoint TSS \
  -b 1000 -a 1000 \
  --binSize 10 \
  -o pergene.all.TSS.matrix.gz \
  



computeMatrix reference-point \
  -S endo_IgG_ratio_w22.bw leaf_IgG_ratio_w22.bw w_endo_H3K27me3_ratio_w22.bw w_endo_H3K56ac_ratio_w22.bw w_leaf_H3K27me3_ratio_w22.bw w_leaf_H3K56ac_ratio_w22.bw \
  -R cluster1_regions.bed cluster2_regions.bed endosperm_teM_genes_regions.bed endosperm_specific_genes_regions.bed pollen_teM_genes_regions.bed core_genes_regions.bed FMEGs_region.bed\
  -p 10 \
  --referencePoint TES \
  -b 1000 -a 1000 \
  --binSize 10 \
  -o pergene.all.TES.matrix.gz 
