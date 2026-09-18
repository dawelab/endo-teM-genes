ml SAMtools/1.21-GCC-13.3.0

module load picard/3.3.0-Java-17


for i in *.sorted.bam; do
base=$(basename "$i" .sorted.bam)
java -jar $EBROOTPICARD/picard.jar MarkDuplicates -I ${i} -O ${base}_nodup.bam -M ${base}_marked_dup_metrics.txt --REMOVE_DUPLICATES true
done

for i in *.sorted.bam ; do
        base=$(basename "$i" .sorted.bam)

        samtools sort -@ 8 -o "${base}.sorted_nodup.bam" "${base}_nodup.bam"
        samtools index -@ 8 "${base}.sorted_nodup.bam"

        samtools view ${base}.sorted_nodup.bam -q 20 -b \
        | samtools sort -@ 8 -o "${base}.sorted_q20.bam" -

  samtools index -@ 8 "${base}.sorted_q20.bam"

done
