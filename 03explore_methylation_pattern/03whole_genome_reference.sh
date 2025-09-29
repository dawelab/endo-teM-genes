#caculate the CG CHG CHH in the whole genome

grep -v "^>" Zm-B73-REFERENCE-NAM-5.0.fa | tr -d '\n' | tr 'a-z' 'A-Z' > genome_clean.txt
fold -w 100 genome_clean.txt > genome_clean_wrapped.txt


# Get total base count
echo -n "Total bases: " > motif_counts.txt
wc -m < genome_clean.txt >> motif_counts.txt

# Count CG
awk 'BEGIN{c=0} {for (i=1; i<=length($0)-1; i++) if (substr($0,i,2)=="CG") c++} END{print "CG:", c}' genome_clean_wrapped.txt >> motif_counts.txt

# Count CHG
awk 'BEGIN{c=0} {
  for (i=1; i<=length($0)-2; i++) {
    if (substr($0,i,1)=="C" && substr($0,i+1,1) ~ /[ATC]/ && substr($0,i+2,1)=="G") c++
  }
} END{print "CHG:", c}' genome_clean_wrapped.txt >> motif_counts.txt

# Count CHH
awk 'BEGIN{c=0} {
  for (i=1; i<=length($0)-2; i++) {
    if (substr($0,i,1)=="C" && substr($0,i+1,1) ~ /[ATC]/ && substr($0,i+2,1) ~ /[ATC]/) c++
  }
} END{print "CHH:", c}' genome_clean_wrapped.txt >> motif_counts.txt

