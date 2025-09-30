# These scripts characterize endosperm teM genes methylation pattern

## 1. Extract the targeted region data into bed file, like 200bp-centered TSS or the gene body regions

## 2. Use bedtools getfasta to get the sequences

## 3. Use CGmapTools to get the mC/total C ratio (C-normalized methylation)

## 4. Count the number of CG, CHG, CHH seperately in the targeted regions

## 5. Calculate the N-normalized methylation based on the algotithm:
$$ N-normalized methylation = \frac{(C-normalized methylation) * (number of cytosine)}{\(total nuclide)} $$