#!/bin/sh

List=$1
Sample_ID=`echo $List | perl -lne 'print $1 if /Abundance\/(.*)/' | sed 's/.\{2\}$//'`
Sample_ID2=`echo $List | perl -lne 'print $1 if /Abundance\/(.*)/' | awk -F '/' '{print $2}'`
gff=$2
wd=$3


awk 'NR==FNR{{a[$1]=$1;next}}{{print $0,a[$1]}}' $List "$gff""$Sample_ID".gff | gawk '$10!=""' | rev | cut -d' ' -f2- | rev > "$List".gff #extracts the lines from the filtered gff that correspond to the list of contigs for each bam
bedtools bamtobed -i "$wd"/"$Sample_ID2".bam > "$wd"/"$Sample_ID2".bam.bed
rm "$wd"/"$Sample_ID2".bam
bedtools coverage -a $List.gff -b "$wd"/"$Sample_ID2".bam.bed -mean > "$List"_coverage #calculates the coverage for each gene in the gff
awk -F '\t|=|;' '{print $1 "_" $10 "\t" $36}' "$List"_coverage #prints the columns of interest


