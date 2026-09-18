ml SRA-Toolkit/3.0.3-gompi-2022a  
ml parallel-fastq-dump/0.6.7-gompi-2022a



for i in {7889760..7889779}; do
    prefetch SRR$i
    fasterq-dump --threads 10 SRR$i -O fastq/
done
