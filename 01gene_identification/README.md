# These scripts identify the different gene sets, including endosperm teM genes, FMEGs, endsoeprm specific genes, and core genes.

## 1. Calclute the TPM from the raw transcript reads  
>Input files:
>1. endo6_38_counts_dec.txt, raw reads generated from the file RNA_Seq_mapping
>2. Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3, B73 genome annotation file
>   
>Algorithm:
>
>>$`\mathrm{TPM}_i = \frac{n_i/l_i}{\sum_{i=1}^j \left(n_i/l_i\right)} \times 10^6`$, $`n_i`$ denotes reads mapped to transcript, $`l_i`$ is the transcript length, j is the total number of genes.  
>>
>Output:
>>all_tpm.txt

## 2. Identify core genes
>Input files:
>
>1. pan_gene_matrix_v3_cyverse.csv downloaded from https://de.cyverse.org/data/ds/iplant/home/shared/NAM/NAM_genome_and_annotation_Jan2021_release/SUPPLEMENTAL_DATA/pangene-files
>2. Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3, B73 genome annotation file
>   
>Criterion:
>>genes are present in all 26 NAM genome lines and the B73 v5 annotation.
>>
>Output:
>>b73_core_genes.txt

## 3.1 identify endosperm teM genes
>Input files:
>1. B73.all.csv downloaded from Yibing Zeng, R Kelly Dawe, Jonathan I Gent, Natural methylation epialleles correlate with gene expression in maize, Genetics(2023).
>2. all_tpm.txt generated from 01calculte_TPM.R

>Criterion:
>1. ≥ 10 CHG and ≥ 40% CG and CHG in the CDS reginos
>2. ≥ 5-fold increase in TPM in endosperm compared to other 9 sporophyte tissues
>3. core genes
>
>Output:
>>endosperm_teM_genes.txt

## 3.2 identify flank-methylated endosperm genes (FMEGs)
>Inpute files:
>1. B73_W22_maternalpreference_with_avg_rpm downladed from Kaitlin Higgins, Vital Nyabashi, Sarah Anderson, Conservation of                                imprinted expression across genotypes is correlated with consistency of imprinting across endosperm development in maize, G3                           Genes|Genomes|Genetics.
 >2. all_tpm.txt generated from 01calculte_TPM.R
 >3. W22_B73_combine.csv, W22 and B73 v5 annotation file 
 >4. Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3, B73 genome annotation file
                    
 >Criterion:
 >1. maternally expressed genes: RPM ≥ 1 & maternal preference ≥ 0.9 
 >2. core genes 
 >3. with DMRs only at the 1kb flank regions
>
>Output:
>>B73_MEGS_FMEGS.txt

## 3.3 identify endosperm specific genes
>Input files:
>>all_tpm.txt generated from 01calculte_TPM.R

>Criterion:
>1. endosperm TPM >=20
>2. ≥ 5-fold increase in TPM in endosperm compared to other 9 sporophyte tissues
>3. core genes
>
>Output:
>>endosperm_specific_genes.txt
                 


     
           
  
