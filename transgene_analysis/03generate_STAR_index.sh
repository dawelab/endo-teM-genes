module load STAR/2.7.10b-GCC-11.3.0

STAR --runThreadN 6 --runMode genomeGenerate --genomeDir /scratch/ys33815/30-1101096464_imprinted-xgene/b73bam/indexgenome --genomeFastaFiles transgene_added_b73_2.fa --sjdbGTFfi
le transgene_added_b73.gff3 --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 100
