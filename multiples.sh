#!/bin/bash
Gene=$1
Matches=$2
Duplicates=$3

if grep -qw ''$Gene'' $Duplicates; then echo -e "$Gene \t Multiple"; else grep -w ''$Gene'' $Matches; fi
