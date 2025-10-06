# Scripts for Characterizing Endosperm teM Gene Methylation Patterns

1. Extract the targeted region data into a BED file, e.g., 200 bp–centered TSS or gene body regions.  
2. Use **bedtools getfasta** to extract the sequences.  
3. Count the number of **CG, CHG, CHH** cytosines separately in the targeted regions.  
4. Use **CGmapTools** to calculate the mC/total C ratio (**C-normalized methylation**).
   Calculate the **N-normalized methylation** using the algorithm:

   $$
   \text{N-normalized methylation} = \frac{\text{C-normalized methylation } \times \text{ number of cytosines}}{\text{total nucleotides}}
   $$

The actrual difference between **C-normalized methylation** and **N-normalized methylation**

    ![two_types_of_methylation_ratio](https://github.com/yiruiS/endo-teM-genes/blob/main/03explore_methylation_pattern/two_types_of_methylation_ratio.png)

5. Calculate the parent specific methylation values.
