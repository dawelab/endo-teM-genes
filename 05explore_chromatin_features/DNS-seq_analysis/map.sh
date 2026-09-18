ml SAMtools/1.21-GCC-13.3.0
ml Bowtie2/2.5.4-GCC-13.3.0

for i in {6898899..6898904} {6898909..6898910}  ; do

  base="SRR${i}"

  bowtie2 --very-sensitive -x b73index -1 "fastq/${base}_1_val_1.fq.gz" -2 "fastq/${base}_2_val_2.fq.gz" -p 16 --rg-id "$base" --rg "SM:$base" --rg "PL:ILLUMINA" \
    | samtools view -@ 8 -bS - \
    | samtools sort -@ 8 -o "${base}.sorted.bam" -

  samtools index -@ 8 "${base}.sorted.bam"

done
