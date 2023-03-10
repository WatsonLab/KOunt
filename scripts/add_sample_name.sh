#!/bin/bash

Fasta=$1
ID=$2

sed -i 's/\s.*$//' $Fasta #remove contig name after first space
sed -i '/>/ s/$/_'$ID'/' $Fasta #add the sample ID to each contig title
