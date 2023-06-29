#!/bin/bash

kegg=$1
cdhit=$2
mmseq=$3

#clusters with mmseq at 90% and 50% PID 
mmseqs createdb "$cdhit""$kegg" "$cdhit""$kegg"_DB    
REASSIGN=TRUE mmseqs cluster "$cdhit""$kegg"_DB "$cdhit""$kegg"_0.9 "$cdhit""$kegg"_tmp_0.9 -v 3 --db-load-mode 3 --cov-mode 3 -c 0.8 --cluster-mode 3  --seq-id-mode 1 --min-seq-id 0.9 --alignment-mode 3 --threads 8 --max-seq-len 60000 --remove-tmp-files 
mmseqs createtsv "$cdhit""$kegg"_DB "$cdhit""$kegg"_DB "$cdhit""$kegg"_0.9 "$mmseq""$kegg"_0.9.tsv
REASSIGN=TRUE mmseqs cluster "$cdhit""$kegg"_DB "$cdhit""$kegg"_0.5 "$cdhit""$kegg"_tmp_0.5 -v 3 --db-load-mode 3 --cov-mode 3 -c 0.8 --cluster-mode 3  --seq-id-mode 1 --min-seq-id 0.5 --alignment-mode 3 --threads 8 --max-seq-len 60000 --remove-tmp-files
mmseqs createtsv "$cdhit""$kegg"_DB "$cdhit""$kegg"_DB "$cdhit""$kegg"_0.5 "$mmseq""$kegg"_0.5.tsv
awk '{print $2 "\t" $1}' "$mmseq""$kegg"_0.9.tsv > "$mmseq""$kegg"_0.9.tsv_tmp
awk 'NR==FNR{a[$1]=$2;next}{print $0,a[$1]}' "$cdhit""$kegg".clstr_reps "$cdhit""$kegg".clstr.txt | awk '{print $3 "\t" $2}' > "$mmseq""$kegg"_0.9_reps
awk 'NR==FNR{a[$1]=$2;next}{print $0,a[$1]}' "$mmseq""$kegg"_0.9.tsv_tmp "$mmseq""$kegg"_0.9_reps | awk '{print $3 "\t" $2}' | sed -e '/^\t*$/d' > "$mmseq""$kegg"_0.9_all
awk '{print $2 "\t" $1}' "$mmseq""$kegg"_0.5.tsv > "$mmseq""$kegg"_0.5.tsv_tmp
awk 'NR==FNR{a[$1]=$2;next}{print $0,a[$1]}' "$cdhit""$kegg".clstr_reps "$cdhit""$kegg".clstr.txt | awk '{print $3 "\t" $2}' > "$mmseq""$kegg"_0.5_reps
awk 'NR==FNR{a[$1]=$2;next}{print $0,a[$1]}' "$mmseq""$kegg"_0.5.tsv_tmp "$mmseq""$kegg"_0.5_reps | awk '{print $3 "\t" $2}' | sed -e '/^\t*$/d' > "$mmseq""$kegg"_0.5_all
mmseqs createsubdb "$cdhit""$kegg"_0.9 "$cdhit""$kegg"_DB "$mmseq""$kegg"_0.9_rep
mmseqs convert2fasta "$mmseq""$kegg"_0.9_rep "$mmseq""$kegg"_0.9.faa
mmseqs createsubdb "$cdhit""$kegg"_0.5 "$cdhit""$kegg"_DB "$mmseq""$kegg"_0.5_rep
mmseqs convert2fasta "$mmseq""$kegg"_0.5_rep "$mmseq""$kegg"_0.5.faa
rm "$mmseq""$kegg"*_0.9_rep*
rm "$mmseq""$kegg"*_0.5_rep*
rm "$mmseq""$kegg"*_0.9.tsv_tmp
rm "$mmseq""$kegg"*_0.5.tsv_tmp
rm "$cdhit""$kegg"*.clstr_reps
rm "$cdhit""$kegg"*.clstr.txt
rm "$cdhit""$kegg"*_DB*
rm "$cdhit""$kegg"_0.9.*
rm "$cdhit""$kegg"_0.5.*
rm -r "$cdhit""$kegg"_tmp_0.9
rm -r "$cdhit""$kegg"_tmp_0.5
