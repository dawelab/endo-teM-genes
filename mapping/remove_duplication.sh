module load picard/2.27.5-Java-15

for f in *Aligned.sortedByCoord.out.bam
do
java -jar $EBROOTPICARD/picard.jar MarkDuplicates -I ${f%%.*}.sortedByCoord.out.bam -O ${f%%.*}_nodup_q20.bam -M marked_dup_metrics.txt --REMOVE_DUPLICATES tr
ue
done
