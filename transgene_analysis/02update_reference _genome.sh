#The transgene.fa file include the transgene sequence
cat Zm-B73-REFERENCE-NAM-5.0.fa transgene.fa > transgene_added_b73.fa
fold -w 80 transgene_added_b73.fa > transgene_added_b73_2.fa

#The transgene.gff3 file treat the transgene as a new chromesome
cat Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 transgene.gff3 > transgene_added_b73.gff3
