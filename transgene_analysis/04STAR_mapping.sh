module load STAR/2.7.10b-GCC-11.3.0

list=("DxJ-25" "DxJ-26" "DxJ-28" "JxD-1" "JxD-2" "JxD-3")

for GENOME in "${list[@]}";
do
STAR --genomeDir ./indexgenome \
--runThreadN 6 \
--readFilesIn ${GENOME}_R1_001_val_1.fq.gz ${GENOME}_R2_001_val_2.fq.gz --readFilesCommand zcat \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFileNamePrefix ./bam/${GENOME}
done
