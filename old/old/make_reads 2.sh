#!/bin/bash
cd "/Users/jdlim/Desktop/RNA/TB/To_individual_genes"
ref=$(<H37Rv.fasta.oneline.txt)
IFS=$';\n'
while read -r gene start end
do
  start="$(echo "$start" | tr -d '[[:space:]]')"
  end="$(echo "$end" | tr -d '[[:space:]]')"
    echo ">${gene}" >> individual_genes_multi.fasta
    cut -c $start-$end H37Rv.fasta.oneline.txt >> individual_genes_multi.fasta
done<RefSeq.cds_H37Rv.csv
