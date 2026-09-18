module load pigz/2.8-GCCcore-13.3.0
ml Python/2.7.18-GCCcore-11.3.0
module load cutadapt/4.9-GCCcore-13.3.0
module load Trim_Galore/0.6.10-GCCcore-12.3.0


for i in {7889827..7889832}; do
  acc="SRR${i}"

trim_galore --fastqc --gzip --paired ${acc}_1.fastq ${acc}_2.fastq -o . --cores 16

done
