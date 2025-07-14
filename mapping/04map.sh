#load packages
module load STAR/2.7.10b-GCC-11.3.0

# go to the file derectory and list all the samples needed
cd /scratch/ys33815/endo6_38/raw
list=("SRR1169626" "SRR1169627" "SRR1169628" "SRR1169629" "SRR1169630" "SRR1169631" "SRR1169632" "SRR1169633" "SRR1169634" "SRR1169635" "SRR1169636" "SRR11696
37" "SRR1169639" "SRR1169640" "SRR1169641" "SRR1169642" "SRR1169643" "SRR1169644" "SRR1169645" "SRR1169646" "SRR1169647")

#mapping
for GENOME in "${list[@]}";
do
STAR --genomeDir /scratch/ys33815/sun/genes/lane2/genomeindex \
--runThreadN 6 \
--readFilesIn ${GENOME}_1_val_1.fq.gz ${GENOME}_2_val_2.fq.gz --readFilesCommand zcat \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFileNamePrefix /scratch/ys33815/endo6_38/bam/${GENOME}
done
