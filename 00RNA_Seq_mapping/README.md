# These scripts mapped the RNA-seq data of maize endosperm from 6 to 38 days after pollination(DAP)
1. RNA-seq data was downloaded from (http://www.ncbi.nlm.nih.gov/sra) under accession number SRP037559.
2. Reads were trimmed using **Trim Galore v0.6.7** with the parameter `--paired`.  
3. The **B73 v5** reference genome was indexed using **STAR**.  
4. Trimmed reads were aligned to the **B73 v5 reference genome** using **STAR**.  
5. Read counts per transcript were quantified using **featureCounts**.  
6. Output the raw read data endo6_38_counts.txt.
