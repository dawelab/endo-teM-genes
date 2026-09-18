ml SAMtools/1.21-GCC-13.3.0
ml Bowtie2/2.5.4-GCC-13.3.0

bowtie2-build Zm-W22-REFERENCE-NRGENE-2.0.fa w22index --threads 10


SRRS=( $(printf "SRR%d\n" {7889828..7889832}) )
base="${SRRS[$SLURM_ARRAY_TASK_ID]}"

TMPDIR="${SLURM_TMPDIR:-/tmp}"
mkdir -p "$TMPDIR/${base}"

bowtie2 --very-sensitive -x w22index \
  -1 "fastq/${base}_1_val_1.fq.gz" \
  -2 "fastq/${base}_2_val_2.fq.gz" \
  -p "$SLURM_CPUS_PER_TASK" \
| samtools sort \
    -@ "$SLURM_CPUS_PER_TASK" \
    -m 2G \
    -T "$TMPDIR/${base}/${base}.tmp" \
    -O BAM \
    -o "${base}.sorted.bam" -

samtools index -@ "$SLURM_CPUS_PER_TASK" "${base}.sorted.bam"
