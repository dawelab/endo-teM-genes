ml SAMtools/1.21-GCC-13.3.0
ml Bowtie2/2.5.4-GCC-13.3.0

bowtie2-build Zm-B73-REFERENCE-NAM-5.0.fa b73index --threads 10

for i in ./fastq/*_trimmed.fq.gz; do

  base=$(basename "$i" _trimmed.fq.gz)

  bowtie2 --very-sensitive -x b73index -U "$i" -p 16 \
    | samtools view -@ 8 -bS - \
    | samtools sort -@ 8 -o "${base}.sorted.bam" -

  samtools index -@ 8 "${base}.sorted.bam"

done
