module load pigz/2.7-GCCcore-11.3.0
ml Python/2.7.18-GCCcore-11.3.0
module load cutadapt/4.5-GCCcore-11.3.0
module load Trim_Galore/0.6.7-GCCcore-11.2.0


list=("DxJ-25" "DxJ-26" "DxJ-28" "JxD-1" "JxD-2" "JxD-3")
for GENOME in "${list[@]}";
do
trim_galore --fastqc --gzip --nextseq-trim 20 --paired ${GENOME}_R1_001.fastq.gz ${GENOME}_R2_001.fastq.gz -o . --cores 4
done
