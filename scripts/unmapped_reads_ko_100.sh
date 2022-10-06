#!/bin/bash
in=$1
wd=$2
best=$3
ref=$4
ko=$5

gawk -F '\t' '$3>=100 {print $1 "\t" $2 "\t" $3 "\t" $8 "\t" $9 "\t" $NF}' "$in" | sort --parallel=8 -k1,1 -k6,6g -k3,3gr | uniq > "$wd"/sorted
awk '! a[$1]++' "$wd"/sorted > "$wd"/sorted_top
awk 'NR==FNR{a[$1];next}$1 in a' <(awk '{print $1 "," $3 "," $6 "\t" $2 "\t" $4 "\t" $5}' "$wd"/sorted_top) <(awk '{print $1 "," $3 "," $6 "\t" $2 "\t" $4 "\t" $5}' "$wd"/sorted) | awk -F ',|\t' 'BEGIN {OFS = "\t"}{print $1,$5,$6,$2,$3,$4}'> "$best" 
rm "$wd"/sorted
rm "$wd"/sorted_top
awk 'BEGIN {OFS = "\t"}{print $1,$2,$3,$6}' "$best" > "$best"_tmp
awk '{print $1}' "$best"_tmp | uniq -c | awk '{print $2 "\t" (1/$1)}' > "$best".bed
awk 'NR==FNR{{a[$1]=$2;next}}{{print $4 "\t" $2 "\t" $3 "\t" a[$1]}}' "$best".bed "$best"_tmp | LC_ALL=C sort --parallel=8 -k1,1 -k2,2n > "$best".bdg
awk 'NR==FNR{{a[$1]=$3;next}}{{print $0,"\t"a[$1]}}' "$ref" "$best".bdg | awk '{print $1 "\t" $2 "\t" $3 "\t" ((($3-$2)/$5)*$4)}' > "$best"_map
awk -F'-' '{print $0 "\t" NF-1}' "$best"_map | awk -v OFMT='%.10f' -v OFS='\t' '{print $1,($4/$5)}' | awk '{gsub(/-/,"\t"$2",")}1' | cut -d, -f2- | sed 's/,$//g' | tr ',' '\n' | awk -v OFMT='%.10f' '{a[$1]+=$2}END{for(i in a) print i,a[i]}' > "$ko"
rm "$best"_tmp
rm "$best".bed
rm "$best".bdg
