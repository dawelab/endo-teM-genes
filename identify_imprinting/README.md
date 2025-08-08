# These scripts describe how do we identify the imprinted genes in maize endoseprm

1. The input file mat_pref_NAM.txt was downloaded from Higgins et al., 2024 "Conservation of imprinted expression across genotypes is correlated with consistency of imprinting across endosperm development in maize"
2. The format transfer file MaizeGDB_maize_pangene_2020_08.tsv was downloaded from https://download.maizegdb.org/Pan-genes/archive/pan-zea-2020/
3. Transfer all the genomes to B73 based on the pangene chart.
4. Select the genes with at least one allele >= 5RPM
5. Group the remaining genes by the B73 gene id and calculte the maternal_preference variance and mean values.
6. Select the genes with <= 0.02 maternal_preference variance, and the maternal_preference means represent the imprinting status.
7. (optional)Draw dot plots showing the maternal_preference distribution of each gene set.
