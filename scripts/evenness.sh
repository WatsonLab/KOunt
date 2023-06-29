#!/bin/sh

List=$1
Sample_ID=`echo $List | perl -lne 'print $1 if /Abundance\/(.*)/' | sed 's/.\{2\}$//'`
Sample_ID2=`echo $List | perl -lne 'print $1 if /Abundance\/(.*)/' | awk -F '/' '{print $2}'`
wd=$2

cut -d';' -f1,14 "$List".gff | cut -f1,9 | sed 's/\tID=/_/g' | sed 's/\;//g' > "$List".contigs #pulls out protein names
bedtools coverage -a "$List".gff -b "$wd"/"$Sample_ID2".bam.bed -d > "$List"_depth #calculates the depth for each gene in the gff
cut -d';' -f1,14 "$List"_depth | cut -f1,9,11 | sed 's/\tID=/_/g' | sed 's/\;//g' > "$List"_depth.2 #pulls out protein id and depth
rm "$List"_depth
groupBy -i "$List"_depth.2 -g 1 -c 2 -o collapse | tr '\t' ',' | datamash transpose -t , --no-strict --filler NA > "$List"_depth.3 #runs bedtools groupby to collapse the depth from each protein
rm "$List"_depth.2
echo "dep <- read.csv(\""$List"_depth.3\", header=TRUE)" > "$List"_depth.R #print header for R script
awk '{print "D=na.omit(dep$"$0"); C=round(mean(D)); D2=D[D<=C]; E=1-(length(D2)-sum(D2)/C)/length(D);E"}' "$List".contigs >> "$List"_depth.R #prints R commands for each protein
R CMD BATCH --vanilla "$List"_depth.R "$List"_depth.Rout #runs R script
rm "$List"_depth.3
grep '^\[1\]' "$List"_depth.Rout | awk '{print $2}' > "$List"_tmp #pulls out evenness
paste "$List".contigs "$List"_tmp > "$List"_evenness #merges protein id and evenness
rm "$List".contigs
rm "$List"_tmp
