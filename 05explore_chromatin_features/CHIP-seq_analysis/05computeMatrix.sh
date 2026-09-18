ml deepTools/3.5.5-gfbf-2023a
ml SAMtools/1.21-GCC-13.3.0

samtools merge -f -o H2AZ_leaf.bam SRR7889760.sorted_q20.bam SRR7889761.sorted_q20.bam 
samtools merge -f -o H3K27ac_leaf.bam SRR7889762.sorted_q20.bam SRR7889763.sorted_q20.bam
samtools merge -f -o H3K27me3_leaf.bam SRR7889764.sorted_q20.bam SRR7889765.sorted_q20.bam
samtools merge -f -o H3K36me3_leaf.bam SRR7889766.sorted_q20.bam SRR7889767.sorted_q20.bam 
samtools merge -f -o H3K4me1_leaf.bam SRR7889768.sorted_q20.bam SRR7889769.sorted_q20.bam
samtools merge -f -o H3K4me3_leaf.bam SRR7889770.sorted_q20.bam SRR7889771.sorted_q20.bam
samtools merge -f -o H3K56ac_leaf.bam SRR7889772.sorted_q20.bam SRR7889773.sorted_q20.bam
samtools merge -f -o H3K9ac_leaf.bam SRR7889774.sorted_q20.bam SRR7889775.sorted_q20.bam
samtools merge -f -o H3_leaf.bam SRR7889776.sorted_q20.bam SRR7889777.sorted_q20.bam
samtools merge -f -o control_leaf.bam SRR7889778.sorted_q20.bam SRR7889779.sorted_q20.bam


for i in *_leaf.bam  ; do
samtools index -@ 10 "$i"
base=$(basename "$i" .bam)
bamCoverage -b ${i} -o ${base}.bw -bs 10 --normalizeUsing RPGC --effectiveGenomeSize 2300000000 -p 10 
done

bigwigCompare \
  -b1 H3_leaf.bw \
  -b2 control_leaf.bw \
  --operation ratio \
  --binSize 10 \
  -o H3_leaf_ratio.bw \
  -p 10

bigwigCompare \
  -b1 H2AZ_leaf.bw \
  -b2 control_leaf.bw \
  --operation ratio \
  --binSize 10 \
  -o H2AZ_leaf_ratio.bw \
  -p 10

bigwigCompare \
  -b1 H3K27ac_leaf.bw \
  -b2 control_leaf.bw \
  --operation ratio \
  --binSize 10 \
  -o H3K27ac_leaf_ratio.bw \
  -p 10

bigwigCompare \
  -b1 H3K27me3_leaf.bw \
  -b2 control_leaf.bw \
  --operation ratio \
  --binSize 10 \
  -o H3K27me3_leaf_ratio.bw \
  -p 10

bigwigCompare \
  -b1 H3K36me3_leaf.bw \
  -b2 control_leaf.bw \
  --operation ratio \
  --binSize 10 \
  -o H3K36me3_leaf_ratio.bw \
  -p 10

bigwigCompare \
  -b1 H3K4me1_leaf.bw \
  -b2 control_leaf.bw \
  --operation ratio \
  --binSize 10 \
  -o H3K4me1_leaf_ratio.bw \
  -p 10

  bigwigCompare \
  -b1 H3K4me3_leaf.bw \
  -b2 control_leaf.bw \
  --operation ratio \
  --binSize 10 \
  -o H3K4me3_leaf_ratio.bw \
  -p 10

  bigwigCompare \
  -b1 H3K56ac_leaf.bw \
  -b2 control_leaf.bw \
  --operation ratio \
  --binSize 10 \
  -o H3K56ac_leaf_ratio.bw \
  -p 10

  bigwigCompare \
  -b1 H3K9ac_leaf.bw \
  -b2 control_leaf.bw \
  --operation ratio \
  --binSize 10 \
  -o H3K9ac_leaf_ratio.bw \
  -p 10

computeMatrix reference-point \
  -S H3_leaf_ratio.bw H2AZ_leaf_ratio.bw H3K27ac_leaf_ratio.bw H3K27me3_leaf_ratio.bw H3K36me3_leaf_ratio.bw H3K4me1_leaf_ratio.bw H3K4me3_leaf_ratio.bw H3K56ac_leaf_ratio.bw H3K9ac_leaf_ratio.bw \
  -R cluster1_regions_b73.bed cluster2_regions_b73.bed endosperm_teM_genes_regions_b73.bed endosperm_specific_genes_regions_b73.bed pollen_teM_genes_regions_b73.bed core_genes_regions_b73.bed FMEGs_region_b73.bed \
  -p 10 \
  --referencePoint TSS \
  -b 1000 -a 1000 \
  --binSize 10 \
  -o pergene.all.TSS.matrix.gz \
 



computeMatrix reference-point \
  -S H3_leaf_ratio.bw H2AZ_leaf_ratio.bw H3K27ac_leaf_ratio.bw H3K27me3_leaf_ratio.bw H3K36me3_leaf_ratio.bw H3K4me1_leaf_ratio.bw H3K4me3_leaf_ratio.bw H3K56ac_leaf_ratio.bw H3K9ac_leaf_ratio.bw \
  -R cluster1_regions_b73.bed cluster2_regions_b73.bed endosperm_teM_genes_regions_b73.bed endosperm_specific_genes_regions_b73.bed pollen_teM_genes_regions_b73.bed core_genes_regions_b73.bed FMEGs_region_b73.bed \
  -p 10 \
  --referencePoint TES \
  -b 1000 -a 1000 \
  --binSize 10 \
  -o pergene.all.TES.matrix.gz 
