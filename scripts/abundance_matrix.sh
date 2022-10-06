#!/bin/bash

Path=$1
KO_List=$2
KEGGS=$3

for file in $Path*KOunt; do cat $file <(awk '{print $1}' $file | cat - $KEGGS | sort | uniq -u | sed 's/$/ 0/') >> "$file"_tmp && mv "$file"_tmp "$file"; done #adds zeros for samples not containing any hits a of particular kegg
for file in $Path*KOunt; do sed -i '1 i\KEGG '${file##*/}'' $file && sed -i 's/_KOunt//' $file; done #adds title line
awk 'BEGIN { OFS = ","} {samples[$1] = samples[$1] OFS $2} END {print "KEGG" samples["KEGG"]; delete samples["KEGG"]; for (KEGG in samples) print KEGG samples[KEGG]}' $Path*KOunt | sed 's/_KOunt//g' > "$Path"hits #merges all of the samples hits
awk '{print $1"\t"substr($0,index($0,$12))}' $KO_List > "$Path"annotation #extracts the kegg id and gene annotation from the ko_list, replacing commas with spaces
awk -F '\t' 'BEGIN { OFS = "\t"} NR==FNR{a[$1]=$2;next}{print $0,a[$1]}' "$Path"annotation $KEGGS | sed 's/,/ /g' | sed 's/\t/,/g' > "$Path"tmp #prints the gene annotation on the keggs list used in matrix generation.
awk -F ',' 'BEGIN { OFS = ","} NR==FNR{a[$1]=$2;next}{print $0,a[$1]}' "$Path"tmp "$Path"hits | awk -F ',' 'BEGIN { OFS = ","}{$2=$NF FS $2;$NF=""}1' | sed 's/\(.*\),/\1/' | sed 's/KEGG,,/KEGG,Description of microbial gene,/g' #adds the annotation to the abundance matrix and moves the last column, gene annotation, to the second column and removes the final comma
rm "$Path"annotation
rm "$Path"hits
rm "$Path"tmp
