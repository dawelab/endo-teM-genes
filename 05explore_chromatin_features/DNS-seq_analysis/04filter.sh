ml SAMtools/1.21-GCC-13.3.0

module load picard/3.3.0-Java-17


bam=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" bam_list.txt)
base=$(basename "$bam" .sorted.RG.bam)
java -jar $EBROOTPICARD/picard.jar MarkDuplicates -I ${bam} -O ${base}_nodup.bam -M ${base}_marked_dup_metrics.txt --REMOVE_DUPLICATES true


samtools sort -@ 8 -o "${base}.sorted_nodup.bam" "${base}_nodup.bam"
samtools index -@ 8 "${base}.sorted_nodup.bam"

samtools view ${base}.sorted_nodup.bam -q 20 -b \
| samtools sort -@ 8 -o "${base}.sorted_q20.bam" -

samtools index -@ 8 "${base}.sorted_q20.bam"
