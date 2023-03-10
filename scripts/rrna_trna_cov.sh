#!/bin/bash

rrna=$1
trna=$2
index=$3

awk 'NR>3{print $1"\t"$3"\t"$4"\t"$5}' "$trna"  > "$trna"_tmp
cat "$trna"_tmp "$rrna" | sort | uniq | awk '{ OFS = "\t"}{if ($2 > $3) print $1,$3,$2,$4; else print $1,$2,$3,$4}' > "$rrna"_cat
awk 'NR==FNR{a[$1]=$2;next}{print $0"\t"a[$4]}' "$index" "$rrna"_cat | gawk '$5>0' | awk '{OFS="\t"}{print $1,$2,$3,$5}' > "$rrna"_catKO
grep 5S "$rrna"_cat | sort -k1,1 -k2,2n > "$rrna"_5S
grep 5.8S "$rrna"_cat | sort -k1,1 -k2,2n > "$rrna"_5.8S
grep 12S "$rrna"_cat | sort -k1,1 -k2,2n > "$rrna"_12S
grep 16S "$rrna"_cat | sort -k1,1 -k2,2n > "$rrna"_16S
grep 18S "$rrna"_cat | sort -k1,1 -k2,2n > "$rrna"_18S
grep 23S "$rrna"_cat | sort -k1,1 -k2,2n > "$rrna"_23S
grep 28S "$rrna"_cat | sort -k1,1 -k2,2n > "$rrna"_28S
grep -w Ala "$rrna"_cat | sort -k1,1 -k2,2n > "$rrna"_Ala
grep -w Arg "$rrna"_cat | sort -k1,1 -k2,2n > "$rrna"_Arg
grep -w Asn "$rrna"_cat | sort -k1,1 -k2,2n > "$rrna"_Asn
grep -w Asp "$rrna"_cat | sort -k1,1 -k2,2n > "$rrna"_Asp
grep -w Cys "$rrna"_cat | sort -k1,1 -k2,2n > "$rrna"_Cys
grep -w Gln "$rrna"_cat | sort -k1,1 -k2,2n > "$rrna"_Gln
grep -w Glu "$rrna"_cat | sort -k1,1 -k2,2n > "$rrna"_Glu
grep -w Gly "$rrna"_cat | sort -k1,1 -k2,2n > "$rrna"_Gly
grep -w His "$rrna"_cat | sort -k1,1 -k2,2n > "$rrna"_His
grep -w Ile "$rrna"_cat | sort -k1,1 -k2,2n > "$rrna"_Ile
grep -w Leu "$rrna"_cat | sort -k1,1 -k2,2n > "$rrna"_Leu
grep -w Lys "$rrna"_cat | sort -k1,1 -k2,2n > "$rrna"_Lys
grep -w Met "$rrna"_cat | sort -k1,1 -k2,2n > "$rrna"_Met
grep -w Phe "$rrna"_cat | sort -k1,1 -k2,2n > "$rrna"_Phe
grep -w Pro "$rrna"_cat | sort -k1,1 -k2,2n > "$rrna"_Pro
grep -w Ser "$rrna"_cat | sort -k1,1 -k2,2n > "$rrna"_Ser
grep -w Thr "$rrna"_cat | sort -k1,1 -k2,2n > "$rrna"_Thr
grep -w Trp "$rrna"_cat | sort -k1,1 -k2,2n > "$rrna"_Trp
grep -w Tyr "$rrna"_cat | sort -k1,1 -k2,2n > "$rrna"_Tyr
grep -w Val "$rrna"_cat | sort -k1,1 -k2,2n > "$rrna"_Val
grep -w SeC "$rrna"_cat | sort -k1,1 -k2,2n > "$rrna"_SeC
grep -w Pyl "$rrna"_cat | sort -k1,1 -k2,2n > "$rrna"_Pyl

if [ -s "$rrna"_cat ]; then bedtools multiinter -i "$rrna"_5S "$rrna"_5.8S "$rrna"_12S "$rrna"_16S "$rrna"_18S "$rrna"_23S "$rrna"_28S "$rrna"_Ala "$rrna"_Arg "$rrna"_Asn "$rrna"_Asp "$rrna"_Cys "$rrna"_Gln "$rrna"_Glu "$rrna"_Gly "$rrna"_His "$rrna"_Ile "$rrna"_Leu "$rrna"_Lys "$rrna"_Met "$rrna"_Phe "$rrna"_Pro "$rrna"_Ser "$rrna"_Thr "$rrna"_Trp "$rrna"_Tyr "$rrna"_Val "$rrna"_SeC "$rrna"_Pyl -names K01985 K01986 K01979 K01977 K01979 K01980 K01982 K14218 K14219 K14220 K14221 K14222 K14223 K14224 K14225 K14226 K14227 K14228 K14229 K14230 K14231 K14232 K14233 K14234 K14235 K14236 K14237 K14238 K14239; else echo "no trna or rrna detected" >&2; fi
rm "$trna"_tmp
rm "$rrna"_Ala
rm "$rrna"_Arg
rm "$rrna"_Asn
rm "$rrna"_Asp
rm "$rrna"_Cys
rm "$rrna"_Gln
rm "$rrna"_Glu
rm "$rrna"_Gly
rm "$rrna"_His
rm "$rrna"_Ile
rm "$rrna"_Leu
rm "$rrna"_Lys
rm "$rrna"_Met
rm "$rrna"_Phe
rm "$rrna"_Pro
rm "$rrna"_Ser
rm "$rrna"_Thr
rm "$rrna"_Trp
rm "$rrna"_Tyr
rm "$rrna"_Val
rm "$rrna"_SeC
rm "$rrna"_Pyl
rm "$rrna"_5S
rm "$rrna"_5.8S
rm "$rrna"_12S
rm "$rrna"_16S
rm "$rrna"_18S
rm "$rrna"_23S
rm "$rrna"_28S
rm "$rrna"_cat
