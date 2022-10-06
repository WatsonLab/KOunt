#!/bin/sh

List=$1
Sample_ID=`echo $List | perl -lne 'print $1 if /Abundance\/(.*)/' | sed 's/.\{2\}$//'`
Sample_ID2=`echo $List | perl -lne 'print $1 if /Abundance\/(.*)/' | awk -F '/' '{print $2}'`
gff=$2
wd=$3

while read file; do awk -v list="$file" '$1 == list' "$gff""$Sample_ID".gff; done < $List > "$List".gff #extracts the lines from the filtered gff that correspond to the list of contigs for each bam
bedtools coverage -a $List.gff -b "$wd"/"$Sample_ID2".bam -mean > "$List"_coverage #calculates the coverage for each gene in the gff
awk -F '\t|=|;' '{print $1 "_" $10 "\t" $36}' "$List"_coverage #prints the columns of interest
