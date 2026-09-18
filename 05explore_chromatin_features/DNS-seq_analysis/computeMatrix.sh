ml deepTools/3.5.5-gfbf-2023a
ml SAMtools/1.21-GCC-13.3.0

samtools merge -f -o EarShoot_HeavyDiges.bam SRR6898899.sorted_q20.bam SRR6898909.sorted_q20.bam -@ 16
samtools merge -f -o EarShoot_LightDiges.bam SRR6898900.sorted_q20.bam SRR6898910.sorted_q20.bam -@ 16
samtools merge -f -o Endo_HeavyDiges.bam SRR6898901.sorted_q20.bam SRR6898903.sorted_q20.bam -@ 16
samtools merge -f -o Endo_LightDiges.bam SRR6898902.sorted_q20.bam SRR6898904.sorted_q20.bam -@ 16

for i in *Diges.bam ; do
  samtools index -@ 16 "$i"
  base=$(basename "$i" .bam)
  bamCoverage -b ${i} -o ${base}.bw -bs 10 --normalizeUsing RPGC --effectiveGenomeSize 2300000000 -p 16 
done

bigwigCompare \
  -b1 EarShoot_LightDiges.bw \
  -b2 EarShoot_HeavyDiges.bw \
  --operation subtract \
  --binSize 10 \
  -o Ear_DNS.bw \
  -p 16

bigwigCompare \
  -b1 Endo_LightDiges.bw \
  -b2 Endo_HeavyDiges.bw \
  --operation subtract \
  --binSize 10 \
  -o Endo_DNS.bw \
  -p 16

computeMatrix reference-point \
  -S Ear_DNS.bw Endo_DNS.bw \
  -R cluster1_regions_b73.bed cluster2_regions_b73.bed endosperm_teM_genes_regions_b73.bed endosperm_specific_genes_regions_b73.bed pollen_teM_genes_regions_b73.bed core_genes_regions_b73.bed FMEGs_region_b73.bed \
  -p 10 \
  --referencePoint TSS \
  -b 1000 -a 1000 \
  --binSize 10 \
  -o pergene.all.TSS.matrix.gz \




computeMatrix reference-point \
  -S Ear_DNS.bw Endo_DNS.bw \
  -R cluster1_regions_b73.bed cluster2_regions_b73.bed endosperm_teM_genes_regions_b73.bed endosperm_specific_genes_regions_b73.bed pollen_teM_genes_regions_b73.bed core_genes_regions_b73.bed FMEGs_region_b73.bed \
  -p 10 \
  --referencePoint TES \
  -b 1000 -a 1000 \
  --binSize 10 \
  -o pergene.all.TES.matrix.gz 
