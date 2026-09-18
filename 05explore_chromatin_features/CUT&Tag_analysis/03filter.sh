ml SAMtools/1.21-GCC-13.3.0

bam=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" bam_list.txt)
base=$(basename "$bam" .sorted.bam)

samtools view ${base}.sorted.bam -q 20 -b \
| samtools sort -@ 8 -o "${base}.q20.sorted.bam" -

samtools index -@ 8 "${base}.q20.sorted.bam"
