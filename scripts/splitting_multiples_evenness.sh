#!/bin/sh

Coverage=$1
KO=$2
WD=$3
Even=$4
PID=$5

#this identifies genes that have hits to >1 KOs and splits the coverage equally between the KOs
#also filters kofamscan results by evenness

awk '{print $1}' $Coverage > "$WD"_1
awk '{print $1}' $KO > "$WD"_2
cat "$WD"_1 "$WD"_2 | sort | uniq -u | awk '{print $1 "\tNoHit"}' | cat "$KO" - > "$WD"_3
awk 'NR==FNR{{a[$1]=$2;next}}{{print $0,a[$1]}}' $Coverage "$WD"_3 > "$WD"_4
cut -f1 "$WD"_4 | sort | uniq -c > "$WD"_5
awk '{print $2 "\t" $1}' "$WD"_5 > "$WD"_6
awk 'NR==FNR{{a[$1]=$2;next}}{{print $0,a[$1]}}' "$WD"_6 "$WD"_4 | awk 'BEGIN { OFS = "\t"} {print $1, $3/$4, $2}' > "$WD"_coverage_tmp
awk 'NR==FNR{{a[$1]=$2;next}}{{print $0,a[$1]}}' "$Even" "$WD"_coverage_tmp | gawk -v pid="$PID" '$4>pid' > "$WD"_coverage
awk '{print $1}' "$WD"_coverage | cat - "$WD"_1 | sort | uniq -u > "$WD"_missing
awk 'NR==FNR{{a[$1]=$2;next}}{{print $0"\t"a[$1]}}' "$WD"_coverage_tmp "$WD"_missing | awk '{print $0"\tNoHit"}' >> "$WD"_coverage
awk '{{print $2 "\t" $3}}' "$WD"_coverage | awk '{{a[$2]+=$1}}END{{for(i in a) print i,a[i]}}' > "$WD"_KOunt
rm "$WD"_1
rm "$WD"_2
rm "$WD"_3
rm "$WD"_4
rm "$WD"_5
rm "$WD"_6
rm "$WD"_missing
rm "$WD"_coverage_tmp
