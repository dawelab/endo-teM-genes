ml deepTools/3.5.5-gfbf-2023a
ml SAMtools/1.21-GCC-13.3.0

samtools merge -f -o atac_leaf.bam SRR7889827.sorted.bam.sorted_q20.bam SRR7889828.sorted.bam.sorted_q20.bam -@ 16

samtools index -@ 16 atac_leaf.bam

bamCoverage -b atac_leaf.bam -o atac_leaf.bw -bs 10 --normalizeUsing RPGC --effectiveGenomeSize 2300000000 -p 16 


bigwigCompare \
  -b atac_leaf.bw \
  --binSize 10 \
  -o atac_leaf.bw \
  -p 16


computeMatrix reference-point \
  -S atac_leaf.bw \
  -R cluster1_regions_b73.bed cluster2_regions_b73.bed endosperm_teM_genes_regions_b73.bed endosperm_specific_genes_regions_b73.bed pollen_teM_genes_regions_b73.bed core_genes_regions_b73.bed FMEGs_region_b73.bed\
  -p 10 \
  --referencePoint TSS \
  -b 1000 -a 1000 \
  --binSize 10 \
  -o pergene.all.TSS.matrix.gz \
  
computeMatrix reference-point \
  -S atac_leaf.bw \
  -R cluster1_regions_b73.bed cluster2_regions_b73.bed endosperm_teM_genes_regions_b73.bed endosperm_specific_genes_regions_b73.bed pollen_teM_genes_regions_b73.bed core_genes_regions_b73.bed FMEGs_region_b73.bed\
  -p 10 \
  --referencePoint TES \
  -b 1000 -a 1000 \
  --binSize 10 \
  -o pergene.all.TES.matrix.gz \
