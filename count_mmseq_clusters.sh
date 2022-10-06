#!/bin/bash

ID=$1

echo "$ID" >> "$ID".txt; awk '{print $1}' "$ID" | sort | uniq -c | awk '{print $1}' | wc -l >> "$ID".txt; awk '{print $1}' "$ID" | sort | uniq -c | awk '{print $1}' | grep -cw 1 >> "$ID".txt; awk '{print $1}' "$ID" | sort | uniq -c | awk '$1>1{c++} END{print c+0}' >> "$ID".txt; wc -l < "$ID" >> "$ID".txt

sed -i -z 's/\n/,/g;s/,$/\n/' "$ID".txt; sed -i 's/_all//g' "$ID".txt
