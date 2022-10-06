#!/bin/bash

kegg_list=$1
fa=$2

perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' "$kegg_list" "$fa" #extracts the genes for each KEGG from the filtered prodigal output

