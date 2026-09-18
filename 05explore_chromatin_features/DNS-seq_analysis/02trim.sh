module load pigz/2.8-GCCcore-13.3.0
ml Python/2.7.18-GCCcore-11.3.0
module load cutadapt/4.9-GCCcore-13.3.0
module load Trim_Galore/0.6.10-GCCcore-12.3.0


for i in {6898899..6898904} {6898909..6898910}; do
trim_galore --fastqc --gzip  --paired SRR${i}_1.fastq SRR${i}_2.fastq -o . --cores 10
done 
