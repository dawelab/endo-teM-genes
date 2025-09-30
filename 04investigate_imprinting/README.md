# Identifying Imprinted Genes in Maize Endosperm

## Input files
> 1. **mat_pref_NAM.txt** — downloaded from Higgins et al., 2024  
>    *"Conservation of imprinted expression across genotypes is correlated with consistency of imprinting across endosperm development in maize"*  
> 2. **MaizeGDB_maize_pangene_2020_08.tsv** — downloaded from [MaizeGDB](https://download.maizegdb.org/Pan-genes/archive/pan-zea-2020/)

## Steps to investigate imprinting status
1. Convert all gene IDs to B73 annotation using `MaizeGDB_maize_pangene_2020_08.tsv`.  
2. Select genes with at least one allele ≥ 5 RPM; discard genes with low expression in both alleles.  
3. Group the remaining genes by B73 gene ID and calculate both the maternal_preference variance and mean.  
4. Select genes with maternal_preference variance ≤ 0.02. The mean maternal_preference value determines the imprinting status.  
5. Visualize results by drawing plots of maternal_preference distribution for each gene set (Fig. 4A, S4A–S4B*).  

