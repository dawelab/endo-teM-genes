module load STAR/2.7.10b-GCC-11.3.0


STAR --runThreadN 6 --runMode genomeGenerate --genomeDir /scratch/ys33815/endo6_38/b73bam/b73index \
--genomeFastaFiles /scratch/ys33815/endo6_38/b73bam/Zm-B73-REFERENCE-NAM-5.0.fa \
--sjdbGTFfile /scratch/ys33815/endo6_38/b73bam/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 100
