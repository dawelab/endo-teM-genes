module load pigz/2.7-GCCcore-11.3.0
ml Python/2.7.18-GCCcore-11.3.0
module load cutadapt/4.5-GCCcore-11.3.0
module load Trim_Galore/0.6.7-GCCcore-11.2.0


list=("SRR1169626" "SRR1169627" "SRR1169628" "SRR1169629" "SRR1169630" "SRR1169631" "SRR1169632" "SRR1169633" "SRR1169634" "SRR1169635" "SRR1169636" "SRR1169637" "SRR1169639" "SRR1169640" "SRR1169641" "SRR1169642" "SRR1169643" "SRR1169644" "SRR1169645" "SRR1169646" "SRR1169647")
for GENOME in "${list[@]}";
do
trim_galore --fastqc --gzip --paired ${GENOME}_1.fastq ${GENOME}_2.fastq -o . --cores 4
done 
