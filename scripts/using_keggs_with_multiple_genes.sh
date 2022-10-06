#!/bin/bash

input=$1
gene=$2

echo $gene | tr ';' '\n' | sed '/^$/d' | awk -v id="$input" '{print id "\t" $0}'
