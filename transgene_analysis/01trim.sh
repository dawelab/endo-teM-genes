module load pigz/2.8-GCCcore-13.3.0
ml Python/3.12.3-GCCcore-13.3.0
module load cutadapt/4.9-GCCcore-13.3.0
module load Trim_Galore/0.6.10-GCCcore-12.3.0


list=("DxJ-25" "DxJ-26" "DxJ-28" "JxD-1" "JxD-2" "JxD-3")
for GENOME in "${list[@]}";
do
trim_galore --fastqc --gzip --nextseq 20 --paired ${GENOME}_R1_001.fastq.gz ${GENOME}_R2_001.fastq.gz -o . --cores 4
done
